#include "OFDM.h"


//////////////////////////////////////
//  FORWARD ERROR CORRECTION (FEC)  //
//////////////////////////////////////

/* Encode 4-bit info block into an 8-bit codeword, by
 * applying a SYSTEMATIC Hamming (8,4) code.
 *
 * Parameters
 *  infoBits - input 4-bit sequence of information bits
 *  encBits  - output 8-bit sequence of encoded bits
 * 
 */
void hamming84_encoder(unsigned char infoBits[], unsigned char encBits[]) {
    encBits[0] = infoBits[0];
    encBits[1] = infoBits[1];
    encBits[2] = infoBits[2];
    encBits[3] = infoBits[3];
    encBits[4] = (infoBits[0]+infoBits[1]+infoBits[3])%2;
    encBits[5] = (infoBits[0]+infoBits[2]+infoBits[3])%2;
    encBits[6] = (infoBits[1]+infoBits[2]+infoBits[3])%2;
    encBits[7] = (infoBits[0]+infoBits[1]+infoBits[2]+infoBits[3]+encBits[4]+encBits[5]+encBits[6])%2;

}


/* Encode input bit sequence by applying Hamming (8,4) code over successive 4-bit blocks.
 *  
 *  dataBits     - input (info) bit sequence
 *  encBits      - output (encoded) bit sequence
 *  infoBitsLen  - input sequence length (num. info bits)
 * 
 */
void fec_encoder(unsigned char infoBits[], unsigned char encBits[], int infoBitsLen) {
    for(int i = 0; i < infoBitsLen/4; i++){
        hamming84_encoder(infoBits + 4*i, encBits + 8*i);
    }
}


////////////////////
//  INTERLEAVING  //
////////////////////

/* Interleave the input sequence by writing cwLen bits into rows,
 * and reading out the obtain matrix column by column.
 *  
 *  bitSeq     - input/output (interleaved) sequence
 *  mat     - interleaving matrix
 *  bitSeqLen     - bit sequence length
 *  numCols - codeword length of the applied FEC code
 * 
 * >>
 * NOTE: Implementation with interleaving done in-place, i.e.
 *       result is placed back to the input array.
 * >>
 */
void bit_interleaver(unsigned char bitSeq[], unsigned char mat[], int bitSeqLen, int numCols) {
    for(int i = 0; i < bitSeqLen; i++){
        mat[i] = bitSeq[i];
    }
    for(int j = 0; j < bitSeqLen; j++){
        int k = (j%(bitSeqLen/numCols))*numCols+floor(j/(bitSeqLen/numCols));
        bitSeq[j] = mat[k];
    }
}

//////////////////////////////////////
//  FORWARD ERROR CORRECTION (FEC)  //
//////////////////////////////////////


/* Decode the SYSTEMATIC Hamming (8,4) code.
 *
 * Parameters
 *  infoBits - input 4-bit sequence of information bits
 *  encBits  - output 8-bit sequence of encoded bits
 * 
 * Output
 *  int - indicator of number of errors, i.e.
 *          0 - no errors
 *          1 - single error (corrected)
 *          2 - double error (uncorrectable)
 * 
 */
int hamming84_decoder(unsigned char encBits[], unsigned char infoBits[]) {
    int s0,s1,s2,s3;
    s0 = (encBits[0]+encBits[1]+encBits[2]+encBits[3]+encBits[4]+encBits[5]+encBits[6]+encBits[7])%2;
    s1 = (encBits[0]+encBits[1]+encBits[3]+encBits[4])%2;
    s2 = (encBits[0]+encBits[2]+encBits[3]+encBits[5])%2;
    s3 = (encBits[1]+encBits[2]+encBits[3]+encBits[6])%2;

    if((s0+s1+s2+s3) == 0)
    {
        infoBits[0] = encBits[0];
        infoBits[1] = encBits[1];
        infoBits[2] = encBits[2];
        infoBits[3] = encBits[3];
        return 0;
    }
    else if ((s1+s2+s3) > 1 && s0 == 1)
    {
        infoBits[0] = encBits[0];
        infoBits[1] = encBits[1];
        infoBits[2] = encBits[2];
        infoBits[3] = encBits[3];
        if(s1+s2+s3 == 3){infoBits[3] = (1+encBits[3])%2;}
        else if(s1+s2 == 2){infoBits[0] = (1+encBits[0])%2;}   
        else if(s1+s3 == 2){infoBits[1] = (1+encBits[1])%2;}
        else if(s2+s3 == 2){infoBits[2] = (1+encBits[2])%2;}
        return 1;
    }
    else{
        infoBits[0] = encBits[0];
        infoBits[1] = encBits[1];
        infoBits[2] = encBits[2];
        infoBits[3] = encBits[3];
        return 2;
    }
    
}


/* Decode bit sequence encoded with Hamming (8,4) code,
 * and report if uncorrectable errors are detected.
 *
 * Parameters
 *  encBits  - input sequence of FEC encoded bits
 *  infoBits - output sequence of (decoded) info bits
 *  infoBitsLen  - output sequence length (num. info bits)
 * 
 * Output
 *  ack_flg - successful transmission acknowledgment
 *              0 - no uncorrectable errors in the input sequence infoBits
 *              1 - at least one codeword is uncorrectable
 */
int fec_decoder(unsigned char encBits[], unsigned char infoBits[], int infoBitsLen) {
    int count = 0, flag = 0;
    for(int i = 0; i < infoBitsLen/4; i++){
        count = hamming84_decoder(encBits + 8*i, infoBits + 4*i);
        if(count == 2){flag = 1;}
    }
    return flag;
}



////////////////////
//  INTERLEAVING  //
////////////////////

/* Deinterleave the input sequence by reverting the interleaving operation.
 *  
 *  bitSeq     - input/output (interleaved) sequence
 *  mat     - deinterleaving matrix
 *  bitSeqLen     - bit sequence length
 *  numCols - number of matrix columns or codeword length of the applied FEC code
 *  
 * >>
 * NOTE: Deinterleaving is done in-place (result is saved to the same array).
 * >>
 */
void bit_deinterleaver(unsigned char bitSeq[], unsigned char mat[], int bitSeqLen, int numCols) {
    for(int i = 0; i < bitSeqLen; i++){
        int k = (i%(bitSeqLen/numCols))*numCols+floor(i/(bitSeqLen/numCols));
        mat[k] = bitSeq[i];
    }
    for(int j = 0; j < bitSeqLen; j++){
        bitSeq[j] = mat[j];
    }
}

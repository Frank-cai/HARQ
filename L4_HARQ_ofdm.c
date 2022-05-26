/*
 This file simulates OFDM with HARQ Type 1 (plain retransmit), where a simple
 Hamming (8,4) code is used for FEC together with interleaving.

 COMPILE with:
    gcc -o L4_HARQ_ofdm L4_HARQ_ofdm.c Transmitter.c Receiver.c Channel.c Auxiliary.c -lm -Wall

 Usage example:
    ./L4_HARQ_ofdm 8 5000

*/

/* Kopic
Compile with command line: 
	make

Usage example:
	./L4_HARQ_ofdm 8 5000
*/

#include "OFDM.h"

void usage(char* progName) {
    printf("\nUsage: %s <taps> <symbol blocks>\n\n", progName);
    printf("taps: length of the channel delay profile\n");
    printf("symbol blocks: number of symbol blocks to simulate\n");
}

void output_file_name(char const* filePath, char const* varname, char* fileName, 
        char* howManyTaps, char* pilotRate) {

    strcpy(fileName, filePath);
    strcat(fileName, varname);
    strcat(fileName, "_HARQ_QAM16_MMSE");

    strcat(fileName, "_L");
    strcat(fileName, howManyTaps);
    strcat(fileName, ".txt");
}


int main(int argc, char** argv) {
    
    // Initialize the random generator. Must be done once per run of main
    srand((unsigned) time(NULL));
    
    // set if SNR-per-bit values for BER/SER evalation
    double EbN0[] = {0, 5, 10, 15, 20, 25, 30, 35, 40};

    // check required if all required parameters are provided
    if (argc != 3) {
        usage(argv[0]);
        return(EXIT_FAILURE);
    }
    
    // set number of bits per symbol (corresponding to QPSK)
    int bitsPerSymbol = QAM64;

    // set equalizer type (MMSE)
    equalizerType equalizer = MMSE;
    
    // read the number of channel taps from parameter list
    int numTaps = atoi(argv[1]);

    // set number of blocks to simulate
    int numBlocks = atoi(argv[2]);


    // FEC parameters > Hammning (8,4) 
    double codeRate = 0.5; // i.e. 4/8
    int cwLen = 8;         // codeword length

    // ARQ parameters
    int max_retrans = 5; // max. number of retransmissions


    // number of bits per block (OFDM symbol)
    int numBits = SYMBOLS_PER_BLOCK * bitsPerSymbol;
    
    // number of information bits (considering the FEC redundancy)
    int numInfoBits = numBits*codeRate;


    // loss due to GI insertion
    double cpLoss = (double) SYMBOLS_PER_BLOCK / (SYMBOLS_PER_BLOCK + CP_LEN);

    // loss due to FEC redundancy
    // NOTE: energy per bit includes both actual bit associated redundancy
    double fecLoss = codeRate;
 
    // Signal magnitude adjustment to match specified EbNo
    // NOTE: Noise is assumed to be of unit power per I/Q component
    double* snr = (double*) malloc(sizeof(EbN0));
    double* sqrtSNR = (double*) malloc(sizeof(EbN0));
    for (int i = 0; i < sizeof(EbN0) / sizeof(double); i++) {
        snr[i] = pow( 10, EbN0[i] / 10); // SNR per bit
        
        // Adjust SNR per symbol to match Eb/No
        snr[i] *= bitsPerSymbol;    // SNR per symbol
        snr[i] *= cpLoss;           // GI (cyclic prefix) insertion
        snr[i] *= fecLoss;          // FEC (redundancy insertion)
        sqrtSNR[i] = sqrt(snr[i]);  // amplitude scaling
    }

    // arrays to keep Tx and Rx information bits in each iteration (i.e. one OFDM symbol)
    unsigned char* txInfoBits = (unsigned char*) malloc( numInfoBits * sizeof(unsigned char));
    unsigned char* rxInfoBits = (unsigned char*) malloc( numInfoBits * sizeof(unsigned char));

    // encoded bits
    unsigned char* txBits = (unsigned char*) malloc( numBits * sizeof(unsigned char));
    unsigned char* rxBits = (unsigned char*) malloc( numBits * sizeof(unsigned char));

    // interleaver and deinterleaver matrix
    unsigned char* intMat = (unsigned char*) malloc( numBits * sizeof(unsigned char));
    unsigned char* deintMat = (unsigned char*) malloc( numBits * sizeof(unsigned char));

    // Tx symbols
    double* txSymI = (double*) malloc( SYMBOLS_PER_BLOCK * sizeof(double));
    double* txSymQ = (double*) malloc( SYMBOLS_PER_BLOCK * sizeof(double));
    
    // Tx OFDM symbol (modulated)
    double* txModI = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    double* txModQ = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    
    // Rx OFDM symbol (after channel)
    double* rxSymI = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    double* rxSymQ = (double*) calloc( SYMBOLS_PER_BLOCK + CP_LEN, sizeof(double));
    
    // Rx symbols
    double* rxEstI = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));
    double* rxEstQ = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));
    
    // multipath channel
    double* hi = (double*) malloc(numTaps * sizeof(double));
    double* hq = (double*) malloc(numTaps * sizeof(double));
    
    double* alpha = (double*) malloc(N * numTaps * sizeof(double));
    double* phi = (double*) malloc(N * numTaps * sizeof(double));

    // Channel transfer function estimate
    double* HestI = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));
    double* HestQ = (double*) malloc((SYMBOLS_PER_BLOCK) * sizeof(double));

    
    // Bit error and retransmission calculation
    double* ber = (double*) malloc(sizeof(EbN0)); 
    double* arq_rate = (double*) malloc(sizeof(EbN0));

    // >> SIMULATION <<

    // ARQ acknowledgment flag (1 errorless frame, 0 erroneous frame)
    unsigned int ackFlg = 1;

    // retransmission counter
    int retrans_cnt = 0;

    // terminal output header
    printf("\nEbNo\tBER\t\tARQ_rate\n");

    // Repeat simulation for each given SNR point
    for (int i = 0; i < sizeof(EbN0) / sizeof(double); i++) {
        
        // initialize bit error and retransmission rate
        ber[i] = 0.0;
        arq_rate[i] = 0.0;

        // Initialize/reset multipath fading simulator
        multipath_fading_init(N, numTaps, fDT, alpha, phi);

        // OFDM simulation loop (block by block)
        for (int j = 0; j < numBlocks; j++) {

            ///////////////////
            //  TRANSMITTER  //
            ///////////////////

            /* YOUR CODE GOES HERE
                       
            // generate information bits or repeat previous frame
            // (in case of requested retransmission)

            // FEC encoding and interleaving

            */
            if (!ackFlg && retrans_cnt < max_retrans){
                retrans_cnt += 1;
            }
            else{
                // generate information bits
                generate_info_bits(numInfoBits, txInfoBits);}

            // FEC encoder
            fec_encoder(txInfoBits, txBits, numInfoBits);
            
            // interleaving
            bit_interleaver(txBits, intMat, numBits, cwLen);

            // modulate bit sequence
            generate_symbols(txBits, bitsPerSymbol, SYMBOLS_PER_BLOCK, txSymI, txSymQ);

            // OFDM modulation
            ifft(txSymI, txSymQ, SYMBOLS_PER_BLOCK, txModI + CP_LEN, txModQ + CP_LEN);
            
            // insert cyclic prefix
            insert_cp(SYMBOLS_PER_BLOCK, CP_LEN, txModI, txModQ);
            

            // scale Tx symbols to match SNR                    
            path_loss(txModI, txModQ, SYMBOLS_PER_BLOCK + CP_LEN, sqrtSNR[i]);
            
            ///////////////
            //  CHANNEL  //
            ///////////////

            // MULTIPATH FADING
            // update channel impulse response at the beginning of each block
            multipath_fading(alpha, phi, N, numTaps, j*(SYMBOLS_PER_BLOCK + CP_LEN), hi, hq);
            
            // Convolution with channel impulse response
            conv_with_gi(hi, hq, numTaps, txModI, txModQ, SYMBOLS_PER_BLOCK + CP_LEN, 
                    CP_LEN, rxSymI, rxSymQ);
            
            // AWGN
            add_noise(rxSymI, rxSymQ, SYMBOLS_PER_BLOCK + CP_LEN);
            
            ////////////////
            //  RECEIVER  //
            ////////////////

            // Automatic Gain Control (AGC)
            // SNR estimation is typically done based on preamble (perfect estimate is assumed here)
            // NOTE: Reverts the pilot scaling introduced in path_loss function
            for (int k = CP_LEN; k < SYMBOLS_PER_BLOCK+CP_LEN; k++) {
                rxSymI[k] /= (sqrtSNR[i]*SQRT_OF_2);
                rxSymQ[k] /= (sqrtSNR[i]*SQRT_OF_2);
            }

            // OFDM demodulation (CI is dropped)
            fft(rxSymI + CP_LEN, rxSymQ + CP_LEN, SYMBOLS_PER_BLOCK, rxEstI, rxEstQ);

            // CTF for assumed perfect channel knowledge
            for (int k = 0; k < SYMBOLS_PER_BLOCK; k++) {
                HestI[k] = 0;
                HestQ[k] = 0;
                for (int m = 0; m < numTaps; m++) {
                    HestI[k] += hi[m] * cos(-2 * PI * k * m / SYMBOLS_PER_BLOCK) 
                            - hq[m] * sin(-2 * PI * k * m / SYMBOLS_PER_BLOCK);
                    HestQ[k] += hq[m] * cos(-2 * PI * k * m / SYMBOLS_PER_BLOCK) 
                            + hi[m] * sin(-2 * PI * k * m / SYMBOLS_PER_BLOCK);
                }
            }

            // Equalization (takes channel transfer function as input)
            fde(equalizer, snr[i], HestI, HestQ, rxEstI, rxEstQ, SYMBOLS_PER_BLOCK);
            
            // Demodulate (detect symbols and output data bits)
            decode_symbols(rxEstI, rxEstQ, SYMBOLS_PER_BLOCK, bitsPerSymbol, rxBits);


            /* YOUR CODE GOES HERE

            // deinterleaving

            // FEC decoding
            
            // request retransmission for unrecoverable errors
            // NOTE: At most max_retrans times for each frame

            */


            // deinterleaving
            bit_deinterleaver(rxBits, deintMat, bitsPerBlock, cwLen);

            // FEC decoder
            if (!fec_decoder(rxBits, rxInfoBits, numBits)){ackFlg = 1;}
            else{ackFlg = 0;}

            if (ackFlg) {
                // count bit and symbol errors in the block
                for (int m = 0; m < numInfoBits; m++) {
                    if (txInfoBits[m] != rxInfoBits[m]) ber[i]++;
                }
            }

        }
            
        // Calculate BER
        ber[i] /= (numInfoBits * numBlocks);

        // Calculate average number of retransmissions
        arq_rate[i] /= (numBlocks);
        
        printf("%.1f\t%f\t%f\n", EbN0[i], ber[i], arq_rate[i]);
    }
    
    // generate output file name for BER
    char fileNameBER[FILE_NAME_SIZE];
    output_file_name("./results/", "BER", fileNameBER, argv[1], argv[3]);
    
    // save BER to file
    save_asignal(EbN0, ber, sizeof(EbN0)/sizeof(double), fileNameBER);
    printf("\nBER saved to %s\n", fileNameBER);
    

    // generate output file name for retranmisssion rate
    char fileNameARQ[FILE_NAME_SIZE];
    output_file_name("./results/", "ARQ", fileNameARQ, argv[1], argv[3]);

    // save retransmission rate to file
    save_asignal(EbN0, arq_rate, sizeof(EbN0)/sizeof(double), fileNameARQ);
    printf("\nARQ rate saved to %s\n\n", fileNameARQ);
    
    
    // Release all allocated memory
    free(txInfoBits);
    free(rxInfoBits);
    free(txBits);
    free(rxBits);
    free(intMat);
    free(deintMat);
    free(txSymI);
    free(txSymQ);
    free(txModI);
    free(txModQ);
    free(rxSymI);
    free(rxSymQ);
    free(rxEstI);
    free(rxEstQ);
    free(snr);
    free(sqrtSNR);
    free(alpha);
    free(phi);
    free(hi);
    free(hq);
    free(ber);
    free(arq_rate);
    //
    free(HestI);
    free(HestQ);
    //
    return EXIT_SUCCESS;
}


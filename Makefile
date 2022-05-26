all:
	gcc -std=c99 -o L4_HARQ_ofdm L4_HARQ_ofdm.c L3_Student.c -I include lib/libOFDM.so -lm -Wall
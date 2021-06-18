/*
 *  example, for timing measurements
 *
 *  MIXMAX, A Pseudo-Random Number Generator
 *
 *  Copyright Konstantin Savvidy.
 *
 *  Free to use, academic or commercial. Do not redistribute without permission.
 *
 *	G.K.Savvidy and N.G.Ter-Arutyunian,
 *  On the Monte Carlo simulation of physical systems,
 *	J.Comput.Phys. 97, 566 (1991);
 *  Preprint EPI-865-16-86, Yerevan, Jan. 1986
 *
 *  K.Savvidy
 *  The MIXMAX random number generator
 *  Comp. Phys. Commun. 196 (2015), pp 161–165
 *  http://dx.doi.org/10.1016/j.cpc.2015.06.003
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "mixmax.h"

// cc -std=c99 -O3 driver_time.c mixmax.c -o time

int main( int argc,  char *argv[] ){
	myuint j,p;
	//double array[ARRAY_SIZE];
	
    fprintf(stderr,"Welcome to the MIXMAX random number generator!\nThe curent matrix size is %u\n"
            "(the actual matrix is not kept in memory in this new efficient implementation)\n"
            "special entry in the matrix is %llu\n"
            "special multiplier is m=2^%u+1\n"
            "Working in the Galois field with modulus 2^%u-1\n", rng_get_N(), (unsigned long long)(SPECIAL), SPECIALMUL, BITS);
	
	fprintf(stderr,"\nHow many numbers?\n(2^10=1024, 2^20=1048576, 2^30=1073741824)\nEnter m: "); 
	if (!scanf("%llu", &p)) {printf("error reading number"); exit(-1);}
	
	rng_state_t S;
	rng_state_t* X = &S;
	// or allocate dynamically:
	//rng_state_t* X=rng_alloc();
	
	//seed_vielbein(X,1);   // seed with unit vector, mostly for testing
	//seed_spbox(X, 123);      // seed with nonlinear SP-box, guarantees a unique seed and independence, but not non-collision of different streams, seed can be 1 ... 2^64-1
	//	seed_lcg(X, 1);        // seed with lcg, guarantees a unique seed and non-collision of different streams, but sadly not independence, seed must be less than 2^61-1
	
	// X->sumtot=apply_bigskip(X->V, X->V, 0, 0, 0, 0); 
	// skip over some steps, here zero -- guarantees complete independence, uniqueness and non-colision of different streams!
	// or just do seed_uniquestream(X, clusterID,  machineID,  runID,   streamID);
    seed_uniquestream(X, 0,  0,  0,  123);

double z=0.0;
	for (j=0; j<p ; j++) { // *ARRAY_SIZE
		z += get_next_float(X); 
		//printf("%1.18F %1.18F %1.18F %1.18F %1.18F\n", get_next_float(X), get_next_float(X), get_next_float(X), get_next_float(X), get_next_float(X) );    
         // for floating point number on [0,1)
		//printf("%1.18F %1.18F %1.18F %1.18F\n", get_next_float(S), get_next_float(S), get_next_float(S), get_next_float(S) );             // for floating point number on [0,1)
		//printf("  %19llu\n", get_next(S) );    // for integer in 0 .. 2^61-2
		//printf("  %19llu\n", j );              // for timing, dont compute anything, just print the counter j
	}
	
	//fprintf(stderr, "%1.18F %1.18F %1.18F %1.18F\n", get_next_float(S), get_next_float(S), get_next_float(S), get_next_float(S) );             // for floating point number on [0,1)
	printf("%1.16F\n", z);
	fprintf(stdout, "ok\n");
	
	// if state was allocated dynamically, free it when its no longer needed:
	//rng_free(S);
	return 0;
}

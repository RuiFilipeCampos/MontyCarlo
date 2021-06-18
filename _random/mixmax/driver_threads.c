/*
 *  example, multithreaded use
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

#include <pthread.h>
#include <stdio.h>
#include "mixmax.h"


#include <stdlib.h>
#include <string.h>
//void *malloc(size_t);
//void *memcpy(void * s1, const void * s2, size_t );
void* go(void *);

unsigned int p;
//rng_state_t* S;

#define MAXTHREADS 32

int main( int argc,  char *argv[] ){
	 int j;
	 int t=1;

	fprintf(stderr,"Welcome to the MIXMAX random number generator!\nThe curent matrix size is %u\n"
			"(the actual matrix is not kept in memory in this new efficient implementation)\n"
			"Working in the Galois field with modulus 2^%u-1\n", rng_get_N(), BITS);
	
	fprintf(stderr,"\nHow many threads [1..%u]?\nEnter t: ", MAXTHREADS); scanf("%d", &t); if (t>MAXTHREADS || t<=0) {fprintf(stderr,"too many threads!\n");return 0;}
	
	fprintf(stderr,"\nHow many numbers?\n(2^10=1024, 2^20=1048576, 2^30=1073741824)\nEnter m: "); scanf("%u", &p);
	
	// allocate dynamically:
	//S=malloc(t*sizeof(rng_state_t));
	pthread_t* thr=malloc(t*sizeof(pthread_t));
	int *tmp=malloc(t*sizeof(int));
	
	setlinebuf(stdout);
	
	for (j=0;j<t;j++){
		int err;
		fprintf(stderr, "Creating thread No %u\n", j);
		 // just passing a variable to the thread does not work : it results in a race condition, instead pass a pointer whose content does not change
		tmp[j]=j;
		err = pthread_create(thr+j, 0,  &go, tmp+j);		
	}
	for (j=0;j<t;j++){
		fprintf(stderr, "Waiting for thread No %u to finish.\n", j);
		pthread_join(*(thr+j), NULL);	
	}
	printf("ok\n");
	
	// if state is allocated dynamically, free it when its no longer needed:
	//rng_free(S);
	free(thr);
	free(tmp);
	return 0;
}

void *go(void *arg){
	int n=*(int *)arg;
	int j;
	char fname[256];
	rng_state_t* X=rng_alloc();
	snprintf(fname, 255, "output_stream_no%03d",n); printf("putting stream No %d into %s\n", n, fname);

	seed_uniquestream(X, 0,0,0, 1+n);
    X->fh=fopen(fname, "w");   //get a file handle for output
    print_state(X); // print the initialized state in the associated file

	for (j=1; j<=p ; j++) {
		fprintf(X->fh, "%1.16F\n", get_next_float(X) );             // for floating point number on [0,1)
	}	
	fprintf(X->fh,"ok\n");
	fclose(X->fh);
	rng_free(X);
	return arg;
}

/*
 *  program for formal verification of correct functioning of the software
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
#include <math.h>
#include "mixmax.h"

#define PRINTONELINE() \
printf("\n{%llu ", X->V[0] ); for (j=1; (j<(Ndim) ); j++) { printf(", %llu ", X->V[j] );} printf("}");

int main( int argc, /* number of arguments */ char *argv[] /* arguments */){
    int j,k;
	int p;

	rng_state_t* X=rng_alloc();
	const int Ndim=rng_get_N();
	
	fprintf(stderr,"Welcome to the MIXMAX random number generator!\nThe curent matrix size is %u\n"
			"(the actual matrix is not kept in memory in this new efficient implementation)\n"
			"special entry in the matrix is %llu\n"
            "special multiplier m=2^%u+1\n"
			"Working in the Galois field with modulus 2^%u-1\n", rng_get_N(), (myuint)SPECIAL, SPECIALMUL, BITS);
    fprintf(stderr,"\nWe will now generate our matrix A, and its power of m\n "); 
    fprintf(stderr,"\nHow many iterations m do you want?\n(2^10=1024, 2^20=1048576, 2^30=1073741824)\n"
			       "Enter m: ");
    
    if(scanf("%d", &p)){fprintf(stderr,"ok\n");}

	printf("\nHere is convenient Mathematica Input:\n AT={");
	for (k=0; (k<Ndim); k++) {
		seed_vielbein(X,k);
		iterate(X);
        PRINTONELINE();
        if (k!=Ndim-1) printf(",");
	}
	printf("} \n A = Transpose[AT] \n");

	printf("\n\n AmT={");
		for (k=0; (k<Ndim); k++) {
			seed_vielbein(X,k);
			for (j=0; (j<p); j++) {
				iterate(X);
			}
            PRINTONELINE();
            if (k!=Ndim-1) printf(",");
		}
		printf("} \n Am = Transpose[AmT] \n"
               // "d = Det[Am] \n cp = CharacteristicPolynomial[AM, x] \n "
			   // "N[Solve[cp == 0, x]] \n N[Eigenvalues[AM]] \n Abs[%%]\n"
				"M1K = Mod[MatrixPower[A, %d], 2^61 - 1]\n"
               "M1K - Am == 0 * IdentityMatrix[%d]\n "
               "Go ahead and run this in Mathematica, M1K should coincide with Am!\n",p,N);
		//\n Maxima input: eq:charpoly(A,x);\n cp:radcan(eq); r:allroots(cp);
	rng_free(X);
	return 0;
}

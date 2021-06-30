/*
 *  example, reading and writing the state to file
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

int main( int argc,  char *argv[] ){
    rng_state_t* X=rng_alloc();
    rng_state_t* Y=rng_alloc();
    
    seed_uniquestream(X, 0,  0,  0,  1);

    FILE* fout;
    if (! (fout = fopen("seed1.conf", "w"))){printf("trouble opening file for writing the state\n");};
    X->fh = fout;
    print_state(X);
    fclose(X->fh );
    read_state(Y, "seed1.conf");
    printf("Successfully read the state from seed1.conf\n");

    int i;
    for(i=0;i<N;i++){
        if (X->V[i] != Y->V[i]) {
            printf("Reading Error: component %d does not match %llu !=%llu\n", i, X->V[i], Y->V[i]);
            exit(-1);
        }
    }
    printf("State Reading Test ok\n");

    rng_free(X);
    rng_free(Y);
	return 0;
}

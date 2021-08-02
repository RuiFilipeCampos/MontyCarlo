/*
 *  driver for testU01 SmallCrush, BigCrush tests 
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

#include "gdef.h"
#include "unif01.h"
#include "sres.h"
#include "bbattery.h"
#include "sknuth.h"
#include "svaria.h"

#include <stdio.h>
#include <stdlib.h>
#include "mixmax.h"

//#define Nstreams 1 // keep it power of two
//rng_state_t X[Nstreams];
rng_state_t S;
rng_state_t* X = &S;

//static int Sindex=0;

/*  // multistream
 uint32_t simpleuint (void){
	Sindex++;
	Sindex &= (Nstreams-1);// Sindex = Sindex % Nstreams;
	return get_next(&X[Sindex]);
} 
*/

//uint64_t A=1;

uint32_t simpleuint (void){
//    return get_next(&X[0]);
    return get_next(X);

} 


double simplefloat (void){return get_next_float(X);} ;

int main (void) {
unif01_Gen *gen;
	setlinebuf(stdout);
	setlinebuf(stderr);

    fprintf(stderr,"Welcome to the MIXMAX random number generator!\nThe curent matrix size is %u\n"
            "(the actual matrix is not kept in memory in this new efficient implementation)\n"
            "special entry in the matrix is %llu\n"
            "special multiplier m=2^%u+1"
            "Working in the Galois field with modulus 2^%u-1\n", rng_get_N(), SPECIAL, SPECIALMUL, BITS);

	//seed_lcg(S, 45634);

	//for (Sindex=0; Sindex < Nstreams; Sindex++){ seed_spbox(&X[Sindex], 1+Sindex);}
    seed_spbox(X, 1);
	
//	gen = unif01_CreateExternGenBits ("uint", simpleuint);
gen = unif01_CreateExternGen01 ("float", simplefloat);
	bbattery_SmallCrush (gen);
//	sknuth_SimpPoker (gen, NULL, 1, 10000000, 0, 64, 64);   // From Crush, MIXMAX for small N<=16 typically fails these
//	sknuth_CouponCollector (gen, NULL, 1, 10000000, 0, 16); // From Crush
//	sknuth_CouponCollector (gen, NULL, 1, 40000000, 0, 4);  // From Crush, N=60 has marginal pvalue=0.998
//	sknuth_Gap(gen, NULL, 1,100000000, 27, 0.0, 0.125);     // Also from Crush
//	//bbattery_Crush (gen); 
	svaria_SumCollector(gen, NULL, 1,500000000, 0, 10.0);      // From BigCrush, some smaller N fail this
	sknuth_Gap(gen, NULL, 1,500000000, 0, 0.0, 0.0625);       // From BigCrush, the only test failed by N=60, N=64 and N=44
	bbattery_BigCrush (gen); 
	unif01_DeleteExternGenBits (gen);

return 0;
}

/*
 * example.cpp
 *
 *  mixmax
 *  A Pseudo-Random Number Generator
 *
 *  Copyright (2008-2016) by Konstantin Savvidy.
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
 *  K.Savvidy and G.Savvidy
 *  Spectrum and Entropy of C-systems. MIXMAX random number generator
 *  Chaos, Solitons & Fractals, Volume 91, (2016) pp. 33–38
 *  http://dx.doi.org/10.1016/j.chaos.2016.05.003
 *
 * 
 */


#include <iostream>
#include <iomanip>
#include <random>

#include "mixmax.hpp"										

int main() {
	mixmax_engine gen{0,0,0,123};		// Create a Mixmax object and initialize the RNG with four 32-bit seeds 0,0,0,1
    //gen.seed(123);  // another way to seed, required by std::random

    fprintf(stderr,"Welcome to the MIXMAX random number generator!\nThe curent matrix size is N=%u\n"
            "(the actual matrix is not kept in memory in this new efficient implementation)\n"
            "special entry in the matrix is %lld\n"
            "special multiplier m=2^%u+1\n"
            "Working in the Galois field with modulus 2^%u-1\n"
            "Generator class size is %lu bytes\n\n",
            gen.rng_get_N(), gen.rng_get_SPECIAL(), gen.rng_get_SPECIALMUL(), 61, sizeof(gen));

      //  gen.print_state();
    
    int n = 10;

//    std::cout << "Print 1200 double precision random numbers:" << std::endl;
//    std::cout << std::fixed; // fixed width
//    for(int i = 0; i < 1200; i++){								// Print 1200 real random numbers
//        std::cout << std::setprecision(18) << std::setw(18)  << gen.get_next_float() << std::endl;
//    }
//    std::cout << "ok" << std::endl;
//
///* Another way to use the generator is through the standard C++11 interface */
//
////    static std::normal_distribution<double> dist{0,1};		// Create a normal distribution
//    static std::uniform_real_distribution<double> dist{0,1};		// Create a uniform distribution
//	
//	std::cout << "Print " << n << " random numbers from given distribution:" << std::endl;
//	
//	for(int i = 0; i < n; i++)								// Print 10 real random numbers
//		std::cout << dist(gen) << std::endl;					

    int p = 1000000000;

    double z=0.0;
    for (int j=0; j<p ; j++) { // *ARRAY_SIZE
        z += gen.get_next_float();
    }
    printf("\n%1.16F\n", z);
    fprintf(stdout, "ok\n");


    return 0;
}

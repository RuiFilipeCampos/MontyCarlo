/*
 Example usage of MIXMAX with GSL
 
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

// cc -std=c99 -O3 -funroll-loops -Wall -I /usr/local/include/ driver_gsl.c mixmax.o -L/usr/local/lib -DHOOKUP_GSL=1  -o gsl -lgsl
// make gsl; GSL_RNG_SEED=123 GSL_RNG_TYPE=MIXMAX ./gsl

#include <stdio.h>
#include "mixmax.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_version.h>

int
main (void)
{
  const gsl_rng_type * T;
  gsl_rng * r;
  unsigned long int seed = 0;

  int i, n = 10;

 // gsl_rng_env_setup();
 // gsl_rng_types_setup ();

  T = gsl_rng_mixmax;
    
    const char *p = getenv ("GSL_RNG_SEED");
    if (p)
    {
        seed = strtoul (p, 0, 0);
        fprintf (stderr, "GSL_RNG_SEED=%lu\n", seed);
    };
    
    gsl_rng_default_seed = seed;

    r = gsl_rng_alloc (T);

    printf("generator type: %s\n", gsl_rng_name (r));
    printf("seed = %lu\n", gsl_rng_default_seed);
    printf("gsl_version=%s\n", gsl_version);
    
  for (i = 0; i < n; i++) 
    {
      double u = gsl_rng_uniform (r);
      printf ("%.16F\n", u);
    }

  gsl_rng_free (r);

  return 0;
}

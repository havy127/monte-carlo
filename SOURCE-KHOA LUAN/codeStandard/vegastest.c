/*
 * This is a sample-program for the vegas-algorithm. Compile the parallel
 * version using something similar to this (once you have the compiled
 * {p|n}vegas object-file):
 *   cc -o vegastest vegastest.c pvegas.o -lm -lpthread  (or equivalent) 
 * and the non-parallel (reference) version with:
 *   cc -o vegastest vegastest.c nvegas.o -lm
 * or the MPI-version equivalently (ask your administrator).  Please make 
 * sure to use the same compiler-switches as in production-runs or the 
 * test will be worthless!
 *
 * The test of the RNG is not really testing the RNG included in the
 * pvegas-sources because they all have different arguments for
 * technical reasons---instead it reemplements it.  It only checks
 * whether they are normal-distributed between 0 and 1.  (It is not the
 * algorithm which is being tested but what the compiler makes out of
 * it.)
 *
 * The test-function is a skewed gaussian that returns 1.0 when
 * integrated over the DIMENSION-dimensional hypercube even if WIDTH
 * gets close to 1.0 because the infinite real axis is mapped on the
 * finite interval (0.0, 1.0).  I do not claim that the choice of
 * vegas-calls is ideal, it is a good start---nothing more. Note
 * also, that the program is not at all suited for measuring timings
 * because testfun() is too lightweight. You need to supply your own
 * functions for that.  Good luck!
 */

#include <math.h>
#include <stdio.h>
#include "vegas.h"
#define DIMENSION 5
#define FUNCTIONS 4
#define RANDTEST 314159
#ifndef PI
#define PI     3.14159265358979323846
#endif

int main(int argc, char **argv)
{

  randomNumberTest();
  

  return(0);
}

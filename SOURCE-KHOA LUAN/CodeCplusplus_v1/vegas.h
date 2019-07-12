
#define SR_P 1279
#define SR_Q 418

#define ALPH 1.5             /* exponent to compress refinement in range     */
#define NDMX 100             /* maximum number of bins per dimension         */
#define MXDIM 6              /* maximum dimension that can be passed by ndim */
#define FNMX 10              /* maximum number of integrands                 */
#define TINY 1.0e-68 /* small, since we are in double-precision      */
#define REPRO 1      /* 0 = default, others used for comparison      */
#define DIMENSION 3
#define FUNCTIONS 4
/* 
 * Different levels of verbosity.  They substitute the old unflexible nprn
 * integer.  Build up your own printlevel by bitwise or-ing several of these.
 */
#define NPRN_INPUT   0x0001  /* print input parameters                       */
#define NPRN_RESULT  0x0002  /* print results of primary integration         */
#define NPRN_SECRES  0x0004  /* print results of secondary integrations      */
#define NPRN_RESULTS 0x0006  /* print results of all integrations            */
#define NPRN_DEFAULT 0x0007  /* corresponds to old default output mode       */
#define NPRN_GRID_8  0x0010  /* print every eighth line of grid data output  */
#define NPRN_GRID_4  0x0020  /* print every fourth line of grid data output  */
#define NPRN_GRID_2  0x0040  /* print every second line of grid data output  */
#define NPRN_GRID    0x0080  /* print all grid data                          */
#define NPRN_ALL     0xffff  /* gossip monger mode: print everything above   */

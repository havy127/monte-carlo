/*
 * NAME
 *   pvegas.c
 *   Parallel version of G.P.Lepage's VEGAS-algorithm.
 *
 * SYNOPSIS
 *   void vegas(double regn[], int ndim, void (*fxn)(double x[], double f[]),
 *              int init, unsigned long ncall, int itmx, int nprn,
 *              int fcns, int pdim, int wrks,
 *              double tgral[], double sd[], double chi2a[]);
 *
 *     regn[]: array specifying the region to be integrated, 2*ndim entries
 *     ndim: dimensionality of space
 *     (*fxn)(x[],f[]): pointer to function to be evaluated (must be MT-safe!)
 *     init: initialization level (start with 0, then 1, later 2)
 *     ncall: number of samples points per iteration
 *     itmx: number of iterations in one call
 *     nprn: bit field, see constants NPRN_* below
 *     fcns: actual number of integrands (<=FNMX),
 *           if > 1 additional function accumulators are used
 *     pdim: dimension of parallel space, 0==autosense, higher=="manual tuning"
 *     wrks: number of parallel working units
 *     tgral[]: pointer to estimate of result (maybe array)
 *     sd[]: pointer to estimate of standard deviation (maybe array)
 *     chi2a[]: pointer to chi-squared over ndim (maybe array)
 *
 *   Compiler-symbols, one of which must be set:
 *     PVEGAS_POSIX     Posix-threads as specified in 1003.1c  (default)
 *     PVEGAS_POSIX_D4  Posix-threads draft 4 (aka old version Pthreads)
 *     PVEGAS_CPS       Convex's CPS-threads
 *
 * DESCRIPTION
 *   pvegas is a parallel farmer-worker implementation of the popular
 *   VEGAS-algorithm. It splits up some dimensions (the first
 *   ndim_par ones) into separate chunks and then lets each
 *   worker-thread evaluate one chunk of these dimensions (and all the
 *   remaining dimensions). The random-numbers from gfsr are
 *   parallelized as well.
 *
 *   This is the multi-threaded version which should compile on most
 *   OSes which support threads (there is another one for MPI-based
 *   systems). It includes Pthreads Draft 4, the final standard (aka
 *   Draft 10) as well as CPS-threads on the Convex Exemplar. It has
 *   been successfully tested on Linux 2.0, Linux 2.2, Linux 2.4,
 *   Solaris 2.5, Solaris 7, DU 4.0, DU 5.0, SPP-UX, OSF/1(3.2),
 *   IRIX 6.4 and IRIX 6.5.  Please feel free to contact
 *   <Richard.Kreckel@Uni-Mainz.DE> if you encounter any
 *   implementation-specific problems.
 *
 *   No external random number generator (RNG) needs to be supplied. A
 *   shift register-generator (SR) is implemented which is initialized
 *   with the system-time. You can force reproducible results by
 *   defining REPRO to be some positive integer different from zero
 *   and sticking to that value. All versions of vegas provided by the
 *   author should return the same numerical values if the same values
 *   for REPRO are used.
 *
 *   Note that the RNG is guaranteed to work properly if your
 *   architecture adheres to the LP64-model, i.e. ints are 32 bit
 *   long. Another, more crucial, assumption is that chars are 8 bit
 *   long. If you are on some strange hardware you are well advised to
 *   check the random numbers manually by consulting the supplied
 *   sample program vegastest.c!
 *
 *   This version may differ considerably from other implementations
 *   found on the net with regards to the arguments passed to it. The
 *   design goal was to have a uniform interface across all versions
 *   of vegas supplied by the author (nvegas.c, pvegas.c,
 *   pvegas_mpi.c) and to make optimal use of parallel facilities.
 *   Consult vegas.h and the samples to make proper use of the
 *   interface.
 *
 * AUTHOR
 * Richard Kreckel, ThEP, Univ. Mainz, October 1996 - October 2000 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "vegas.h"

#if defined(PVEGAS_CPS)         /* Definitions for CPS-threads: */
#define MUTEX_T cps_mutex_t
#define MUTEX_LOCK cps_mutex_lock
#define MUTEX_UNLOCK cps_mutex_unlock
#define MUTEX_DESTROY cps_mutex_free
#define THREAD_EXIT cps_thread_exit()
#include <cps.h>
#elif defined(PVEGAS_POSIX_D4)  /* Definitions for old-style Posix-threads: */
#define MUTEXATTR_DEF pthread_mutexattr_default
#define THREADATTR_DEF pthread_attr_default
#define MUTEX_T pthread_mutex_t
#define MUTEX_LOCK pthread_mutex_lock
#define MUTEX_UNLOCK pthread_mutex_unlock
#define MUTEX_DESTROY pthread_mutex_destroy
#define THREAD_EXIT pthread_exit(dummy)
#include <pthread.h>
#else                           /* Definitions for Pthreads (POSIX 1003.1c) */
#define MUTEXATTR_DEF NULL      /* (The default) */
#define THREADATTR_DEF NULL
#define MUTEX_T pthread_mutex_t
#define MUTEX_LOCK pthread_mutex_lock
#define MUTEX_UNLOCK pthread_mutex_unlock
#define MUTEX_DESTROY pthread_mutex_destroy
#define THREAD_EXIT pthread_exit(dummy)
#include <pthread.h>
#endif

#define TINY 1.0e-68         /* small, since we are in double-precision      */
#define MXWORK 32            /* no more than MXWORK parallel threads, please */
#define REPRO 1              /* 0 = default, others used for comparison      */

MUTEX_T ndim_par_mtx;        /* lock to prevent double eval. of slices       */
MUTEX_T rewrite_mtx;         /* lock to keep rewrite of results ordered      */
static int mds;              /* ==1: statified smpl. (see README)  (mds)     */
int nd;                      /* slices in grid (c.f. NDMX)         (nd)      */
int ng;                      /*                                    (ng)      */
int npg;                     /* number of calls within bin         (npg)     */
int ndim_par;                /* dimensionality of parallel space             */
int gndim;                   /* global copy of ndim                          */
double dxg;                  /*                                    (dxg)     */
double xnd;                  /*                                    (xnd)     */
double xJac;                 /* Jacobian of integration            (xjac)    */
typedef struct {
  double ti;                 /* sum for f over bins                (ti)      */
  double tsi;                /* sum for variances over bins        (tsi)     */
} binAccu;                   /* accumulator over bins / hypercubes...        */
binAccu Ab[FNMX];            /* ...one for each integrand                    */
double d[NDMX][MXDIM];       /*                                    (d[][])   */
double di[NDMX][MXDIM];      /* delta i                            (di[][])  */
double dx[MXDIM];            /* width of integration region        (dx[])    */
double gregn[MXDIM];         /* global copy of regn                          */
static double xi[MXDIM][NDMX]; /*                                  (xi[][])  */
int kgl[MXDIM-1];            /* first few dimensions are stored here         */
void (*p_fxn)(double [], double []);  /* copy of fxn                         */
unsigned int gfsr_m[MXWORK][SR_P];  /* status-field for GFSR-generators      */
int gfsr_k[MXWORK];          /* pointer into status-fields                   */
int gfsr_n[MXWORK];          /* number of independent generator              */
static int gfsr_not_initialized = 1;  /* flag: need to initialize GFSR-field */
unsigned int rdum;           /* linear congruential counter in kw_rand()     */
double gfsr_norm;            /* will be set such that gfsr is normalized     */
int functions;               /* copy of (*ctl).fcns                          */

/*
 * This routine is used by vegas. It rebins a vector of densities xi into
 * new bins defined by a vector r.
 */
void rebin(double rc, int nd, double r[], double xin[], double xi[])
{
  int i;
  int k = 0;
  double dr = 0.0;
  double xn = 0.0;
  double xo = 0.0;
  
  for (i=0; i<nd-1; i++) {
    while (rc > dr) {
      dr += r[k++];
    }
    if (k > 1) xo = xi[k-2];
    xn = xi[k-1];
    dr -= rc;
    xin[i] = xn-(xn-xo)*dr/r[k-1];
  }
  for (i=0; i<nd-1; i++) xi[i] = xin[i];
  xi[nd-1] = 1.0;
}

/*
 * gfsr produces the random numbers once the starting values have been
 * initialized using gfsr_init.
 */
double gfsr_rand(unsigned int w[], int *k)
{
  int j;
  
  (*k)++;
  if (*k >= SR_P) *k = 0;
  j = *k + SR_Q;
  if (j >= SR_P) j -= SR_P;
  w[*k] = w[*k] ^ w[j];
  return((double)w[*k] * gfsr_norm);
}

/*
 * Simple linear congruential generator used to initialize SR.
 * The multiplier and increment were provided by KOMA, Univ.Mainz.
 * (We may be abusing them a little, but I have checked that 
 * it has maximal periodicity which is what we want.)
 */
unsigned int kw_rand(void)
{
  rdum = (1812433253*rdum + 314159265);
  return (rdum);
}

/*
 * gfsr_init initializes the sequences using values from kw_rand.
 */
void gfsr_init(long seed)
{
  int i, j;
  
  printf("Initializing SR-sequences with seed %ld\n",seed);
  gfsr_norm = (double)(1/(pow(2.0,(long double)(8*sizeof(int))) - 1));
  rdum = (unsigned int)seed;
  for (i=0; i<MXWORK; i++) {
#if (REPRO != 0)
    rdum = (unsigned int)seed;
#endif
    for (j=0; j<SR_P; j++)
      gfsr_m[i][j] = kw_rand();  /* initialize starting values               */
    gfsr_k[i] = -1;          /* initialize pointer                           */
    gfsr_n[i] = i;           /* set current number                           */
  }
  gfsr_not_initialized = 0;
}

/*
 * This function rewrites the results collected by one worker into
 * the global variables. It must be locked (using rewrite_mtx)!
 */
void p_vegasrewrite(binAccu r_Ab[FNMX],
                    double r_d[NDMX][MXDIM], double r_di[NDMX][MXDIM])
{
  int i, j;
  
  for (j=0; j<functions; j++) {
    Ab[j].ti += r_Ab[j].ti;
    Ab[j].tsi += r_Ab[j].tsi;
  }
  for (j=0; j<gndim; j++) {
    for (i=0; i<nd; i++) {
      d[i][j] += r_d[i][j];
      di[i][j] += r_di[i][j];
    }
  }
}

/*
 * The implementation of a VEGAS-worker with its local variables
 * and its code.
 */
void* p_vegasloop(void *dummy)
{
  int i, j, k, h, l;
  int ia[MXDIM];             /*                                    (ia[])    */
  int kg[MXDIM];             /*                                    (kg[])    */
  typedef struct {
    double f2;               /* f squared                          (f2)      */
    double fb;               /* sum for f within bin               (fb)      */
    double f2b;              /* sum for f2 within bin              (f2b)     */
    unsigned long npg;       /* number of calls within bin f != 0            */
  } pointAccu;               /* accumulator over points x within bins...     */
  pointAccu Ax[FNMX];        /* ...one for each integrand                    */
  binAccu p_Ab[FNMX];        /* accumulator for bins (parallel space only!)  */
  double f[FNMX];            /* array passed into fxn for evaluation at x    */
  double x[MXDIM];           /* evaluation point                   (x[])     */
  double rc;                 /*                                    (rc)      */
  double wgt;                /* weight                             (wgt)     */
  double p_xnd;              /* local copy of xnd                            */
  double xn;                 /*                                    (xn)      */
  double xo;                 /*                                    (xo)      */
  double p_xJac;             /* local copy of xJac                           */
  double p_d[NDMX][MXDIM];   /* local copy of d[][]                          */
  double p_di[NDMX][MXDIM];  /* local copy of di[][]                         */
  int p_ng;                  /* local copy of ng                             */
  int p_npg;                 /* local copy of npg                            */
  int p_ndim;                /* local copy of ndim                           */
  double p_ndim_par;         /* local copy of ndim_par                       */
  double p_mds;              /* local copy of switch mds -- don't touch      */
  double p_dxg;              /* local copy of dxg                            */
  double p_dx[MXDIM];        /* local copy of dx[]                           */
  double p_xi[MXDIM][NDMX];  /* local copy of xi[][]                         */
  double p_regn[MXDIM];      /* local copy of regn                           */
  int* gfsr_arg = (int *)dummy;  /* this number will be _our_ generator      */
#if (REPRO != 0)
  int kgl_rep[MXDIM-1];      /* copy of kgl, used only if if REPRO != 0      */
  unsigned long i_rep;       /* number of RNs to skip                        */
  int f_rep = 0;             /* a flag: have we iterated yet?                */
#endif
  
  p_ndim = gndim;
  p_ndim_par = ndim_par; 
  p_ng = ng;
  p_npg = npg; 
  p_mds = mds;
  p_dxg = dxg; 
  p_xnd = xnd; 
  p_xJac = xJac;
  for (i=p_ndim_par; i<p_ndim; i++) kg[i] = 1;
  if (p_ndim_par==0 && *(int *)dummy>0) THREAD_EXIT;  /* makes things foolproof */
  for (j=0; j<p_ndim; j++) {
    p_dx[j] = dx[j];
    p_regn[j] = gregn[j];
    for (i=0; i<nd; i++) p_xi[j][i] = xi[j][i];
  }
  
  for (j=0; j<functions; j++) {
    p_Ab[j].ti = p_Ab[j].tsi = 0.0;
  }
  for (j=0; j<p_ndim; j++) {
    kg[j] = 1;
    for (i=0; i<nd; i++) p_d[i][j] = p_di[i][j] = 0.0;
  }
  
#if (REPRO != 0)
  for (i=0; i<p_ndim_par; i++) kgl_rep[i] = 1;
#endif
  for (;;) {  /* Here the outer iteration begins, the one which is split up */
    if (MUTEX_LOCK(&ndim_par_mtx)!=0)
      perror("Mutex ndim_par_mtx on fire");
    if (kgl[0]==0 &&  /* we are done already! */
        p_ndim_par>0) /* may not make much sense, but makes things foolproof */
      {
        if (MUTEX_UNLOCK(&ndim_par_mtx)!=0)
          perror("Mutex ndim_par_mtx on fire");
        if (MUTEX_LOCK(&rewrite_mtx)!=0)
          perror("Mutex rewrite_mtx on fire");
        p_vegasrewrite(p_Ab,p_d,p_di);
        if (MUTEX_UNLOCK(&rewrite_mtx)!=0)
          perror("Mutex rewrite_mtx on fire");
#if (REPRO != 0)
        if (f_rep) i_rep = 0; else {
          for (k=1,i=0; i<p_ndim-p_ndim_par; i++) k *= p_ng;
          i_rep = k*npg*p_ndim;
        }
        for (j=0; j<p_ndim_par; j++) {
          for (k=1,i=0; i<p_ndim-j-1; i++) k *= p_ng;
          i = (p_ng-kgl_rep[j]);
          i_rep += k*i*npg*p_ndim;
        }
        while (i_rep-- > 0)
          gfsr_rand(gfsr_m[*gfsr_arg],&gfsr_k[*gfsr_arg]);
#endif
        THREAD_EXIT;
      }
    for (i=0; i<p_ndim_par; i++) kg[i] = kgl[i];
    for (l=p_ndim_par; l>=1; l--) {
      kgl[l-1] %= p_ng;
      if (++kgl[l-1] != 1) break;
    }
    if (l < 1) kgl[0] = 0;  /* we are done! */
#if (REPRO != 0)
    i_rep = 0;
    for (j=0; j<p_ndim_par; j++) {
      for (k=1,i=0; i<p_ndim-j-1; i++) k *= p_ng;
      i = (kg[j]-kgl_rep[j]);
      if (j==p_ndim_par-1 && f_rep) i--;
      i_rep += k*i*npg*p_ndim;
    }
    while (i_rep-- > 0) 
      gfsr_rand(gfsr_m[*gfsr_arg],&gfsr_k[*gfsr_arg]);
    
    for (i=0; i<p_ndim_par; i++) kgl_rep[i] = kg[i];
#endif
    if (MUTEX_UNLOCK(&ndim_par_mtx)!=0)
      perror("Mutex ndim_par_mtx on fire");
    for (;;) {  /* ...and here the inner one, which is untouched */
      for (j=0; j<functions; j++) {
        Ax[j].fb = 0.0;
        Ax[j].f2b = 0.0;
        Ax[j].npg = 0;
      }
      for (k=0; k<p_npg; k++) {
        wgt = p_xJac;
        for (j=0; j<p_ndim; j++) {
          xn = (kg[j]-gfsr_rand(gfsr_m[*gfsr_arg],&gfsr_k[*gfsr_arg]))*p_dxg+1.0;
          ia[j] = ((int)xn<NDMX) ? ((int)xn) : (NDMX);
          ia[j] = (ia[j]>1) ? (ia[j]) : (1);
          if (ia[j] > 1) {
            xo = p_xi[j][ia[j]-1]-p_xi[j][ia[j]-2];
            rc = p_xi[j][ia[j]-2]+(xn-ia[j])*xo;
          } else {
            xo = p_xi[j][ia[j]-1];
            rc = (xn-ia[j])*xo;
          }
          x[j] = p_regn[j]+rc*p_dx[j];
          wgt *= xo*p_xnd;
        }
        (*p_fxn)(x,f);  /* call integrand at point x */
        for (j=0; j<functions; j++) {
          if (f[j] != 0.0) ++Ax[j].npg;
          f[j] *= wgt;
          Ax[j].f2 = f[j]*f[j];
          Ax[j].fb += f[j];
          Ax[j].f2b += Ax[j].f2;
        }
        for (j=0; j<p_ndim; j++) {
          p_di[ia[j]-1][j] += f[0];
          if (p_mds >= 0) p_d[ia[j]-1][j] += Ax[0].f2;
        }
      }  /* end of loop within hypercube */
      for (j=0; j<functions; j++) {
        Ax[j].f2b = sqrt(Ax[j].f2b*Ax[j].npg);
        Ax[j].f2b = (Ax[j].f2b-Ax[j].fb)*(Ax[j].f2b+Ax[j].fb);
        if (Ax[j].f2b <= 0.0) Ax[j].f2b = TINY;
        p_Ab[j].ti += Ax[j].fb;
        p_Ab[j].tsi += Ax[j].f2b;
      }
      if (p_mds < 0) {
        for (j=0; j<p_ndim; j++) p_d[ia[j]-1][j] += Ax[0].f2b;
      }
      for (h=p_ndim-1; h>=p_ndim_par; h--) {
        kg[h] %= p_ng;
        if (++kg[h] != 1) break;
      }
      if (h < p_ndim_par) break;
    }  /* end of loop over parallel-space hypercubes */
#if (REPRO != 0)
    f_rep = 1;
#endif
    if (l < 1) break;
  }  /* end of loop over orthogonal-space hypercubes */
  if (MUTEX_LOCK(&rewrite_mtx)!=0)
    perror("Mutex rewrite_mtx on fire");
  p_vegasrewrite(p_Ab,p_d,p_di);
  if (MUTEX_UNLOCK(&rewrite_mtx)!=0)
    perror("Mutex rewrite_mtx on fire");
  THREAD_EXIT;
  return NULL; /* stops complaints of IRIX-compiler */
}

/*
 * The routine pvegas to be called by the user. Parameters are just like
 * in vegas, except that workers specifies the degree of parallelization.
 */
void vegas(double regn[], int ndim, void (*fxn)(double x[], double f[]),
           int init, unsigned long ncall, int itmx, int nprn,
           int fcns, int pdim, int wrks,
           double tgral[], double sd[], double chi2a[])
{
#ifdef PVEGAS_CPS
  int worker[MXWORK];
  int which, node;
#else
  pthread_t worker[MXWORK];
#endif
  static int ndo;            /*                                    (ndo)     */
  int it;                    /* iteration counter                  (it)      */
  static int ittot;          /* iteration counter across init>1              */
  int i, j, k;               /* counters                           (i, j, k) */
  double calls;              /* real total number of calls to fxn  (calls)   */
  double dv2g;               /*                                    (dv2g)    */
  double rc;                 /*                                    (rc)      */
  double xn;                 /*                                    (xn)      */
  double xo;                 /*                                    (xo)      */
  double xin[NDMX];          /* aux. variable for rebinning        (xin[])   */
  typedef struct {
    double Wgt;              /* weight                             (wgt)     */
    double sWgt;             /* cumulative sum for weights         (swgt)    */
    double sChi;             /* cumulative sum for chi^2           (schi)    */
    double sInt;             /* cumulative sum for integral        (si)      */
  } iterAccu;                /* accumulator for stuff at end of iteration... */
  static iterAccu Ai[FNMX];  /* ...one for each integrand                    */
  double dt[MXDIM];          /*                                    (dt[])    */
  double r[NDMX];            /*                                    (r[])     */
  int wrks2;                 /* counts workers spawned successfully          */
  int wmax;                  /* limit, if workers too large                  */
  
  wrks = (wrks<1) ? 1 : ((MXWORK<wrks) ? MXWORK : wrks);
  gndim = ndim;
  if (pdim == 0)
    ndim_par = ndim/2;
  else
    ndim_par = (pdim<ndim) ? (pdim) : (ndim-1);
  functions = (fcns<FNMX) ? (fcns) : (FNMX);
  
  for (i=0; i<ndim_par; i++) kgl[i] = 1;
  p_fxn = fxn;
  for (j=0; j<ndim; j++) gregn[j] = regn[j];
#ifdef PVEGAS_CPS
  if (cps_mutex_alloc(&rewrite_mtx)!=0)
    perror("Mutex rewrite_mtx drowned");
  if (cps_mutex_alloc(&ndim_par_mtx)!=0)
    perror("Mutex ndim_par_mtx drowned");
#else
  if (pthread_mutex_init(&rewrite_mtx,MUTEXATTR_DEF)!=0)
    perror("Mutex rewrite_mtx drowned");
  if (pthread_mutex_init(&ndim_par_mtx,MUTEXATTR_DEF)!=0)
    perror("Mutex ndim_par_mtx drowned");
#endif
  
#if (REPRO != 0)
  if (gfsr_not_initialized) gfsr_init(REPRO);
#else
  if (gfsr_not_initialized) gfsr_init((long)time(NULL));
#endif
  if (init <= 0) {    /* entry for cold start        */
    mds = ndo = 1;    /* Careful with that Axe, Eugene!    *
                       * mds=0 will trash parallelization. */
    for (j=0; j<ndim; j++) xi[j][0] = 1.0;
  }
  if (init <= 1) {    /* inherit the previous grid   */
    for (j=0; j<functions; j++) {
      Ai[j].sInt = 0.0;
      Ai[j].sWgt = 0.0;
      Ai[j].sChi = 0.0;
    }
    ittot = 1;
  }
  if (init <= 2) {    /* inherit grid and results    */
    nd = NDMX;
    ng = 1;
    if (mds) {
      ng = (int)pow(ncall/2.0+0.25,1.0/ndim);
      mds = 1;
      if ((2*ng-NDMX) >= 0) {
        mds = -1;
        npg = ng/NDMX+1;
        nd = ng/npg;
        ng = npg*nd;
      }
    }
    wmax = 1; 
    for (i=0; i<ndim_par; i++) wmax *= ng;
    if (wrks>wmax) wrks = wmax;
    for (k=1,i=0; i<ndim; i++) k *= ng;
    npg = (ncall/k>2) ? (ncall/k) : (2);
    calls = (double)npg * (double)k;
    dxg = 1.0/ng;
    for (dv2g=1,i=0; i<ndim; i++) dv2g *= dxg;
    dv2g = calls*calls*dv2g*dv2g/npg/npg/(npg-1.0);
    xnd = nd;
    dxg *= xnd;
    xJac = 1.0/calls;
    for (j=0; j<ndim; j++) {
      dx[j] = regn[j+ndim]-regn[j];
      xJac *= dx[j];
    }
    if (nd != ndo) {
      for (i=0; i<(nd>ndo?nd:ndo); i++) r[i] = 1.0;
      for (j=0; j<ndim; j++) rebin(ndo/xnd,nd,r,xin,xi[j]);
      ndo = nd;
    }
    if (nprn & NPRN_INPUT) {
      printf("%s:  ndim= %3d  ncall= %8.0f     %3d thread(s)\n",
             " Input parameters for vegas",ndim,calls,wrks);
      printf("%28s  ittot=%5d  itmx=%5d    %5d^%1d hypercubes\n"," ",ittot,itmx,ng,ndim_par);
      printf("%28s  nprn=0x%04x  ALPH=%5.2f\n"," ",nprn,ALPH);
      printf("%28s  mds=%3d  nd=%4d%15s npg=%d\n"," ",mds,nd," ",npg);
      for (j=0; j<ndim; j++) {
        printf("%30s xl[%2d]= %11.4g xu[%2d]= %11.4g\n",
               " ",j,regn[j],j,regn[j+ndim]);
      }
    }
  }
  
  for (it=ittot; it<=itmx+ittot-1; it++) {
    wrks2 = wrks;
    kgl[0] = 1;  /* kgl[0] is also used to determine if we are done yet */
    for (i=0; i<ndim_par; i++) 
      if (kgl[i]!=1) {
        fprintf(stderr,"Warning: kgl[%d] not correctly initialized!\n",i);
        kgl[i] = 1;
      }
    for (i=0; i<functions; i++)
      Ab[i].ti = Ab[i].tsi = 0.0;
    for (j=0; j<gndim; j++) {
      for (i=0; i<nd; i++) d[i][j] = di[i][j] = 0.0;
    }
#if defined (PVEGAS_CPS)
    node = CPS_ANY_NODE;
    for(i=0; i<wrks; i++) {
      if ((k = cps_thread_create(&node,(void (*)(void *))(p_vegasloop),(void *)&gfsr_n[i])<0)) {
        fprintf(stderr,"cps_thread_create=%d\n",k);
        perror("Thread has drowned");
        wrks2--;
      }
    }
    if (wrks2 == 0) {
      fprintf(stderr,"Uargh: All our workers drowned at creation!\n");
      exit(-1);
    }
    if (wrks2 != wrks) 
      fprintf(stderr,"Only %d out of %d workers could be created.\n",wrks2,wrks);
    which = 1;
    for(i=0; i<wrks2; i++) {
      if (cps_thread_wait(&which)==-1)
        perror("Thread went Rambo");
    }
#else
    for(i=0; i<wrks; i++) {
      if (pthread_create(&worker[i],THREADATTR_DEF,
                         p_vegasloop,
                         (void *)&gfsr_n[i]) != 0) {
        perror("Thread has drowned");
        wrks2--;
      }
    }
    if (wrks2 == 0) {
      fprintf(stderr,"Uargh: All our workers drowned at creation!\n");
      exit(-1);
    }
    if (wrks2 != wrks)
      fprintf(stderr,"Only %d out of %d workers could be created.\n",wrks2,wrks);
    for(i=0; i<wrks2; i++) {
      if (pthread_join(worker[i],NULL)!=0)
        perror("Thread went Rambo");
    }
#endif
#if (REPRO != 0)                   /* Now advance the remaining RNGs in case */
    for (j=wrks2; j<MXWORK; j++) {  /* the next call runs on more threads */
      for (i=0; i<SR_P; i++)
        gfsr_m[j][i] = gfsr_m[0][i];
      gfsr_k[j] = gfsr_k[0];
    }
#endif
    for (j=0; j<functions; j++) {
      Ab[j].tsi *= dv2g;
      Ai[j].Wgt = 1.0/Ab[j].tsi;
      Ai[j].sInt += Ai[j].Wgt*Ab[j].ti;
      Ai[j].sChi += Ai[j].Wgt*Ab[j].ti*Ab[j].ti;
      Ai[j].sWgt += Ai[j].Wgt;
      tgral[j] = Ai[j].sInt/Ai[j].sWgt;
      chi2a[j] = (Ai[j].sChi-Ai[j].sInt*tgral[j])/(it-0.9999);
      if (chi2a[j] < 0.0) chi2a[j] = 0.0;
      sd[j] = sqrt(1.0/Ai[j].sWgt);
      Ab[j].tsi = sqrt(Ab[j].tsi);
    }
    if (nprn & NPRN_RESULT) {
      printf("%s %3d : integral = %14.7g +/-  %9.2g\n",
             " iteration no.",it,Ab[0].ti,Ab[0].tsi);
      printf("%s integral =%14.7g+/-%9.2g  chi^2/IT n = %9.2g\n",
             " all iterations:  ",tgral[0],sd[0],chi2a[0]);
    }
    if (nprn & NPRN_SECRES) {
      for (i=1; i<functions; i++) {
        printf("   %4d%s%14.7g+/-%9.2g  chi^2/IT n = %9.2g\n",
               i,".additional integral= ",tgral[i],sd[i],chi2a[i]);
      }
    }
    if (nprn & (NPRN_GRID | NPRN_GRID_2 | NPRN_GRID_4 | NPRN_GRID_8)) {
      for (j=0; j<ndim; j++) {
        printf(" data for axis  %2d\n",j);
        printf("%6s%13s%11s%13s%11s%13s\n", 
               "X","delta i","X","delta i","X","delta i");
        for (i=0; i<nd; i += 3) {
          for (k=0; k<3 && i+k<nd; k++) {
            printf("%8.5f%12.4g    ",xi[j][i+k],di[i+k][j]);
          }
          printf("\n");
          if (nprn & NPRN_GRID_8) k = 3*(8-1);
          if (nprn & NPRN_GRID_4) k = 3*(4-1);
          if (nprn & NPRN_GRID_2) k = 3*(2-1);
          if (nprn & NPRN_GRID) k = 3*(1-1);
          i += k;
        }
      }
    }
    if (nprn) fflush(NULL);
    for (j=0; j<ndim; j++) {
      xo = d[0][j];
      xn = d[1][j];
      d[0][j] = (xo+xn)/2.0;
      dt[j] = d[0][j];
      for (i=1; i<nd-1; i++) {
        rc = xo+xn;
        xo = xn;
        xn = d[i+1][j];
        d[i][j] = (rc+xn)/3.0;
        dt[j] += d[i][j];
      }
      d[nd-1][j] = (xo+xn)/2.0;
      dt[j] += d[nd-1][j];
    }
    for (j=0; j<ndim; j++) {
      rc = 0.0;
      for (i=0; i<nd; i++) {
        if (d[i][j] < TINY) d[i][j] = TINY;
        r[i] = pow((1.0-d[i][j]/dt[j])/
                   (log(dt[j])-log(d[i][j])),ALPH);
        rc += r[i];
      }
      rebin(rc/xnd,nd,r,xin,xi[j]);
    }
  }
  ittot += itmx;
  if (MUTEX_DESTROY(&rewrite_mtx)!=0)
    perror("Mutex rewrite_mtx went Rambo");
  if (MUTEX_DESTROY(&ndim_par_mtx)!=0)
    perror("Mutex ndim_par_mtx went Rambo");
}

#undef MUTEX_T 
#undef MUTEX_LOCK
#undef MUTEX_UNLOCK
#undef MUTEX_DESTROY
#undef SR_P
#undef SR_Q
#undef REPRO
#undef MXWORK
#undef TINY

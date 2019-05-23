/*
 * NAME
 *   eduvegas.c
 *   Implementation of G.P.Lepage's VEGAS-algorithm and
 *      Richard Kreckel, ThEP, Univ. Mainz, October 1996 - October 2000 
 *
 *
 * AIM eduvegas.c
 *   modified for education purposes: Explain statistics knowlegde in vegas
 *
 * STATUS
 *   Under modifing .... current time: November 2017,
 *
 * SYNOPSIS
 *   void vegas(double regn[], int ndim, void (*fxn)(double x[], double f),
 *              int init, unsigned long ncall, int itmx, int nprn,
 *              double tgral[], double sd[], double chi2a[]);
 *
 *     regn[]: array specifying the region to be integrated, 2*ndim entries
 *     ndim: dimensionality of space
 *     (*fxn)(x[],f[]): pointer to function to be evaluated (must be MT-safe!)
 *     init: initialization level (start with 0, then 1, later 2)
 *     : number of samples points per iteration
 *     itmx: number of iterations in one call
 *     nprn: bit field, see constants NPRN_* below
 *     tgral[]: pointer to estimate of result (maybe array)
 *     sd[]: pointer to estimate of standard deviation (maybe array)
 *     chi2a[]: pointer to chi-squared over ndim (maybe array)*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "vegas.h"

#define TINY 1.0e-68 /* small, since we are in double-precision      */
#define REPRO 1      /* 0 = default, others used for comparison      */

unsigned int gfsr_m[SR_P];           /* status-field for the GFSR-generator          */
int gfsr_k;                          /* pointer into status-field                    */
static int gfsr_not_initialized = 1; /* flag: need to initialize GFSR-field */
unsigned int rdum;                   /* linear congruential counter in kw_rand()     */
double gfsr_normn;                   /* will be set such that gfsr is normalized     */
int functions;                       /* copy of (*ctl).fcns                          */

/*
 * This routine is used by vegas. It rebins a vector of densities xi into
 * new bins defined by a vector r. Origin from vegas
 */
/*
Điều chỉnh lại các độ dài delta xi 
*/
void rebin(double wgt_avg, int nd, double wgt[], double xin[], double xi[])
{
  int i;
  int k = 0;
  double residues = 0.0;
  double xn = 0.0;
  double xo = 0.0;

  for (i = 0; i < nd - 1; i++)
  {
    while (wgt_avg > residues)
    {
      residues += wgt[k++];
    }
    if (k > 1)
      xo = xi[k - 2];
    xn = xi[k - 1];
    residues -= wgt_avg; //residues of accummulate weights up to bin k-1
    xin[i] = xn - (xn - xo) * residues / wgt[k - 1];
  }
  for (i = 0; i < nd - 1; i++)
    xi[i] = xin[i];
  xi[nd - 1] = 1.0;
}

/*
 * gfsr produces the random numbers once the starting values have been
 * initialized using gfsr_init. Origin from vegas
 */

/*
* Hàm tạo mãng ngẫu nhiên. 
*/
double gfsr_rand(unsigned int w[], int *k)
{
  int j;

  (*k)++;
  if (*k >= SR_P)
    *k = 0;
  j = *k + SR_Q;
  if (j >= SR_P)
    j -= SR_P;
  w[*k] = w[*k] ^ w[j];
  return ((double)w[*k] * gfsr_normn);
}

/**  Origin from vegas
 * Simple linear congruential generator used to initialize SR.
 * The multiplier and increment were provided by KOMA, Univ.Mainz.
 * (kw stands for Kalos & Whitlock)
 */
/*
* Giã số ngẫu nhiên
*/
unsigned int kw_rand(void)
{
  rdum = (1812433253 * rdum + 314159265);
  return (rdum);
}

/**  Origin from vegas
 * gfsr_init initializes the sequences using values from kw_rand.
 */
void gfsr_init(long seed)
{
  int i, j;

  printf("Initializing SR-sequences with seed %ld\n", seed);
  gfsr_normn = (double)(1 / (pow(2.0, (long double)(8 * sizeof(int))) - 1));
  rdum = (unsigned int)seed;
  for (j = 0; j < SR_P; j++)
    gfsr_m[j] = kw_rand(); /* initialize starting values */
  gfsr_k = -1;             /*         initialize pointer */
  gfsr_not_initialized = 0;
}
//=======================================================

void randomNumberTest()
{
  
   double num;
    gfsr_init((long)time(NULL));
    num = gfsr_rand(gfsr_m, &gfsr_k);
   printf ("Random number: %f",num); 

}
//=======================================================
void vegas(double regn[], int ndim, void (*fxn)(double x[], double *f), //replace f[] by f(1)
           int init, unsigned long ncall, int itmx, int nprn,           // fxn(x, &f); /* call integrand at point x */
           double *tgralin, double *sdin, double *chi2ain)
{
  double tgral, sd, chi2a;
  static int ndo;   /*                                    (ndo)     */
  int it;           /* iteration counter                  (it)      */
  static int ittot; /* iteration counter across init>1              */
  int i, j, k;      /* counters                           (i, j, k) */
  int nd;           /* slices in grid (c.f. NDMX)         (nd)      */
  int ng;           /* real number of bins on an axis     (ng)      */
  unsigned int npg; /* number of calls within bin         (npg)     */
  static int mds;   /* ==1: statified smpl.               (mds)     */
  int ia[MXDIM];    /* index of bin (ng bin per axis)     (ia[])    */
  int kg[MXDIM];    /* contral npg calls locate same bin  (kg[])    */
  double calls;     /* real total number of calls to fxn  (calls)   */
  double dv2g;      /*                                    (dv2g)    */
  double dxg;       /* ratio of number of grids per bins  (dxg)     */
  double rc;        /*                                    (rc)      */
  double wgt;       /* weight                             (wgt)     */
  double xn;        /*                                    (xn)      */
  double xnd;       /* another name of nd                 (xnd)     */
  double xo;        /*                                    (xo)      */
  double xJac;      /* Jacobian of integration            (xjac)    */
  typedef struct
  {
    double Wgt;       /* weight                             (wgt)     */
    double sWgt;      /* cumulative sum for weights         (swgt)    */
    double sChi;      /* cumulative sum for chi^2           (schi)    */
    double sInt;      /* cumulative sum for integral        (si)      */
  } iterAccu;         /* accumulator for stuff at end of iteration... */
  static iterAccu Ai; /* ...for integrand                             */
  typedef struct
  {
    double ti;  /* sum for f over bins                (ti)      */
    double tsi; /* sum for variances over bins        (tsi)     */
  } binAccu;    /* accumulator over bins / hypercubes...        */
  binAccu Ab;   /* ...for integrand                             */
  typedef struct
  {
    double f2;                   /* f squared                          (f2)      */
    double fb;                   /* sum for f within bin               (fb)      */
    double f2b;                  /* sum for f2 within bin              (f2b)     */
    unsigned long npg;           /* number of calls within bin f != 0            */
  } pointAccu;                   /* accumulator over points x within bins...     */
  pointAccu Ax;                  /* ...for integrand                             */
  double f;                      /* passed into fxn for evaluation at x          */
  double x[MXDIM];               /* evaluation point                   (x[])     */
  double d[NDMX][MXDIM];         /*                                    (d[][])   */
  double di[NDMX][MXDIM];        /* delta i                            (di[][])  */
  double dt[MXDIM];              /*                                    (dt[])    */
  double r[NDMX];                /*                                    (r[])     */
  static double xi[MXDIM][NDMX]; /*                                    (xi[][])  */
  double xin[NDMX];              /* aux. variable for rebinning        (xin[])   */
  double dx[MXDIM];              /* width of integration region        (dx[])    */
  double xrand;                  /* uniform dens 0.0 <= xrand < 1.0   (xrand)   */

#if (REPRO != 0)
  if (gfsr_not_initialized)
    gfsr_init(REPRO);
#else
  if (gfsr_not_initialized)
    gfsr_init((long)time(NULL));
#endif
  if (init <= 0)
  { /* entry for cold start        */
    mds = ndo = 1;
    for (j = 0; j < ndim; j++)
      xi[j][0] = 1.0;
  }
  if (init <= 1)
  { /* inherit the previous grid   */
    Ai.sInt = 0.0;
    Ai.sWgt = 0.0;
    Ai.sChi = 0.0;
    ittot = 1;
  }
  if (init <= 2)
  { /* inherit grid and results    */
    nd = NDMX;
    ng = 1;
    if (mds)
    {
      ng = (int)pow(ncall / 2.0 + 0.25, 1.0 / ndim); // M=2N^n; ng <-- N
      mds = 1;                                       //concentrate subintervals where integrand is largest in magnitude
      if ((2 * ng - NDMX) >= 0)
      {
        mds = -1; //concentrate subintervals where integrand is largest in error (sigma)
        npg = ng / NDMX + 1;
        nd = ng / npg;
        ng = npg * nd;
      }
    }
    for (k = 1, i = 0; i < ndim; i++)
      k *= ng; // k <-- ng^ndim
    npg = (ncall / k > 2) ? (ncall / k) : (2);
    calls = (double)npg * (double)k;
    dxg = 1.0 / ng;
    for (dv2g = 1, i = 0; i < ndim; i++)
      dv2g *= dxg;
    dv2g = calls * calls * dv2g * dv2g / npg / npg / (npg - 1.0); //dv2g <-- (npg.ng^ndim)^2/ngp^2/ng^(2ndim)/(npg-1)
    xnd = nd;
    dxg *= xnd; // dxg <-- nd/ng
    xJac = 1.0 / calls;
    for (j = 0; j < ndim; j++)
    {
      dx[j] = regn[j + ndim] - regn[j];
      xJac *= dx[j];
    }
    if (nd != ndo)
    { //ndo=1 at this step: initialize bin
      for (i = 0; i < (nd > ndo ? nd : ndo); i++)
        r[i] = 1.0;
      for (j = 0; j < ndim; j++)
        rebin(ndo / xnd, nd, r, xin, xi[j]);
      ndo = nd; //number of subdivisions on an axis
    }
    if (nprn & NPRN_INPUT)
    {
      printf("%s:  ndim= %3d  ncall= %8.0f\n",
             " Input parameters for vegas", ndim, calls);
      printf("%28s  ittot=%5d  itmx=%5d\n", " ", ittot, itmx);
      printf("%28s  nprn=0x%04x  ALPH=%5.2f\n", " ", nprn, ALPH);
      printf("%28s  mds=%3d  nd=%4d%15s npg=%d\n", " ", mds, nd, " ", npg);
      for (j = 0; j < ndim; j++)
      {
        printf("%30s xl[%2d]= %11.4g xu[%2d]= %11.4g\n",
               " ", j, regn[j], j, regn[j + ndim]);
      }
    }
  }
  for (it = ittot; it <= itmx + ittot - 1; it++)
  {
    Ab.ti = Ab.tsi = 0.0;
    for (j = 0; j < ndim; j++)
    {
      kg[j] = 1;
      for (i = 0; i < nd; i++)
        d[i][j] = di[i][j] = 0.0;
    }
    for (;;)
    {
      Ax.fb = 0.0; //refresh to start a hypercube
      Ax.f2b = 0.0;
      Ax.npg = 0;
      for (k = 0; k < npg; k++)
      {
        wgt = xJac;
        for (j = 0; j < ndim; j++)
        {
          xrand = gfsr_rand(gfsr_m, &gfsr_k); //uniform random number, prob <--1/ng
          xn = (kg[j] - xrand) * dxg + 1.0;   //1<kg[j]<=ng
          ia[j] = ((int)xn < NDMX) ? ((int)xn) : (NDMX);
          ia[j] = (ia[j] > 1) ? (ia[j]) : (1);
          if (ia[j] > 1)
          {
            xo = xi[j][ia[j] - 1] - xi[j][ia[j] - 2];
            rc = xi[j][ia[j] - 2] + (xn - ia[j]) * xo;
          }
          else
          {
            xo = xi[j][ia[j] - 1];
            rc = (xn - ia[j]) * xo;
          }
          x[j] = regn[j] + rc * dx[j];
          wgt *= xo * xnd; // ??? p(xrand)dxrand = p(rc)drc  ???, where p(xrand) = 1/ng
        }
        fxn(x, &f); /* call integrand at point x */
        if (f != 0.0)
          ++Ax.npg;
        f *= wgt; // contribution of x to I
        Ax.f2 = f * f;
        Ax.fb += f; //sum over npg points whithin hypercube-->contribution of an hypercube to I
        Ax.f2b += Ax.f2;
        for (j = 0; j < ndim; j++)
        {
          di[ia[j] - 1][j] += f;
          if (mds >= 0)
            d[ia[j] - 1][j] += Ax.f2; //sum up npg magnitude f2 whithin each bin of a hypercube
        }
      } /* end of loop within hypercube */
      Ax.f2b = sqrt(Ax.f2b * Ax.npg);
      Ax.f2b = (Ax.f2b - Ax.fb) * (Ax.f2b + Ax.fb); //sum{aCube}{f^2}-(sum{aCube}f)^2-->sigma^2 acube
      if (Ax.f2b <= 0.0)
        Ax.f2b = TINY;
      Ab.ti += Ax.fb;   //sum over hypercubes --> sum{over calls points}f --> I of an iteration
      Ab.tsi += Ax.f2b; //sum(cubes){sigma^2 of cube) --> sigtot^2 of an iteration
      if (mds < 0)
      {
        for (j = 0; j < ndim; j++)
          d[ia[j] - 1][j] += Ax.f2b; //sum(over cubs per bin){sigma^2)-->sigbin^2
      }
      for (k = ndim - 1; k >= 0; k--)
      {
        kg[k] %= ng; //create kg[k]: 1<=kg[k]<=ng; ng^ndim values (equal or none)
        if (++kg[k] != 1)
          break;
      }
      if (k < 0)
        break;
    } /* end of loop over hypercubes */
    Ab.tsi *= dv2g;
    Ai.Wgt = 1.0 / Ab.tsi;             // ~ 1/(sigtot_i)^2 of iteration ith
    Ai.sInt += Ai.Wgt * Ab.ti;         // ~ sum(iteration){I_i/(sigtot_i)^2}
    Ai.sChi += Ai.Wgt * Ab.ti * Ab.ti; // ~ sum(iteration){(I_i)^2/(sigtot_i)^2}
    Ai.sWgt += Ai.Wgt;                 // ~ sum(iteration){1/(sigtot_i)^2}-->1/sigI^2
    tgral = Ai.sInt / Ai.sWgt;         //~sigI^2.sum(iteration){I_i/(sigtot_i)^2}
    chi2a = (Ai.sChi - Ai.sInt * tgral) / (it - 0.9999);
    if (chi2a < 0.0)
      chi2a = 0.0;
    sd = sqrt(1.0 / Ai.sWgt); //sigI^2
    Ab.tsi = sqrt(Ab.tsi);
    *sdin = sd;
    *tgralin = tgral;
    *chi2ain = chi2a;
    if (nprn & NPRN_RESULT)
    {
      printf("%s %3d : integral = %14.7g +/-  %9.2g\n",
             " iteration no.", it, Ab.ti, Ab.tsi);
      printf("%s integral =%14.7g+/-%9.2g  chi^2/IT n = %9.2g\n",
             " all iterations:  ", tgral, sd, chi2a);
    }
    if (nprn & (NPRN_GRID | NPRN_GRID_2 | NPRN_GRID_4 | NPRN_GRID_8))
    {
      for (j = 0; j < ndim; j++)
      {
        printf(" data for axis  %2d\n", j);
        printf("%6s%13s%11s%13s%11s%13s\n",
               "X", "delta i", "X", "delta i", "X", "delta i");
        for (i = 0; i < nd; i += 3)
        {
          for (k = 0; k < 3 && i + k < nd; k++)
          {
            printf("%8.5f%12.4g    ", xi[j][i + k], di[i + k][j]);
          }
          printf("\n");
          if (nprn & NPRN_GRID_8)
            k = 3 * (8 - 1);
          if (nprn & NPRN_GRID_4)
            k = 3 * (4 - 1);
          if (nprn & NPRN_GRID_2)
            k = 3 * (2 - 1);
          if (nprn & NPRN_GRID)
            k = 3 * (1 - 1);
          i += k;
        }
      }
    }
    if (nprn)
      fflush(NULL);
    for (j = 0; j < ndim; j++)
    {
      xo = d[0][j];              //sigb0
      xn = d[1][j];              //sigb1
      d[0][j] = (xo + xn) / 2.0; //(sigb0+sigb1)/2
      dt[j] = d[0][j];           //(sigb0+sigb1)/2
      for (i = 1; i < nd - 1; i++)
      {
        rc = xo + xn;
        xo = xn;
        xn = d[i + 1][j];
        d[i][j] = (rc + xn) / 3.0; //(sigb{i-1}+sigbi+sigb{i+1})/3
        dt[j] += d[i][j];
      }
      d[nd - 1][j] = (xo + xn) / 2.0; //(sigb{nd-1}+sigbnd)/2
      dt[j] += d[nd - 1][j];          //(sigb0+sigb1)/2+..(sigb{i-1}+sigbi+sigb{i+1})/3+..(sigb{nd-1}+sigbnd)/2
    }
    for (j = 0; j < ndim; j++)
    {
      rc = 0.0;
      for (i = 0; i < nd; i++)
      {
        if (d[i][j] < TINY)
          d[i][j] = TINY;
        r[i] = pow((1.0 - d[i][j] / dt[j]) /
                       (log(dt[j]) - log(d[i][j])),
                   ALPH); //sigbi contribution of bin ith in sumsigbin
        rc += r[i];
      }
      rebin(rc / xnd, nd, r, xin, xi[j]); //rebin wrt input ratio sigbi/sumsigbin
    }
  }
  ittot += itmx;
}

#undef SR_P
#undef SR_Q
#undef REPRO
#undef TINY

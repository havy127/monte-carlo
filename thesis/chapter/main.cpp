
#define ALPH 1.5             /* exponent to compress refinement in range     */
#define NDMX 100             /* maximum number of bins per dimension         */
#define MXDIM 6              /* maximum dimension that can be passed by ndim */
#define TINY 1.0e-68 /* small, since we are in double-precision      */
#define DIMENSION 1

#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include <stdio.h>
#include "vegas.h"
using namespace std;


class info_iteration
{
    public:
    double Wgt;     /* wgt = 1/(sigma)^2 của mỗi vòng lặp           */
    double sWgt;    /* Tổng wgt                                     */
    double sInt;    /* Tổng ước tính tích phân qua từng vòng lặp    */
};

class info_bins 
{
    public:
    double ti;  /* tổng ước tính tích phân của tổng các ô trong 1 vòng lặp        */
    double tsi; /* Tổng phương sai của tổng các ô trong 1 vòng lặp                */
    
};

class info_bin
{   
    public:
    double f2;                 /* f^2 -> f bình phương                                 */
    double fb;                 /* Tổng f bên trong 1 bin. Mỗi bin được gieo 2 điểm     */
    double f2b;                /* Tổng f2 bên trong 1 bin.                             */
    unsigned long samples_2;   /* số điểm gieo bên trong 1 bin                         */
    
};

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
        residues -= wgt_avg; 
        xin[i] = xn - (xn - xo) * residues / wgt[k - 1];
    }
    for (i = 0; i < nd - 1; i++)
    {
        xi[i] = xin[i];
    }
    xi[nd - 1] = 1.0;
}

double gfsr_rand (double a, double b){

   return a + (b - a)*rand()/RAND_MAX;

}

void vegas(double regn[], int ndim, void (*fxn)(double x[], double *f), 
           int init, unsigned long samples, int itmx,          
           double *tgralin, double *sdin)
{
    double tgral, sd;
    static int ndo;  
    int it;           
    static int ittot; 
    int i, j, k;      
    int nd;           
    int ng;           
    unsigned int samples_2; 
    static int mds;  
    int index[MXDIM];   
    int kg[MXDIM];    
    double calls;     
    double dv2g;     
    double dxg;       
    double rc,rc1;        
    double wgt;       
    double xn,xn1;    
    double xnd;       
    double xo,xo1;       
    double xJac;      
    static info_iteration Ai;
    static info_bins Ab;    
    static info_bin Ax;                
    double f;                     
    double x[MXDIM];              
    double d[NDMX][MXDIM];         
    double dt[MXDIM];             
    double r[NDMX];               
    static double xi[MXDIM][NDMX];
    double xin[NDMX];              
    double dx[MXDIM];              
    double xrand;                  

    if (init <= 0)
    { 
        mds = ndo = 1;
        for (j = 0; j < ndim; j++)
            xi[j][0] = 1.0;
    }
    if (init <= 1)
    { 
        Ai.sInt = 0.0;
        Ai.sWgt = 0.0;
        ittot = 1;
    }
    if (init <= 2)
    {
        nd = NDMX;
        ng = 1;
        if (mds)
        {
            ng = (int)pow(samples/ 2.0 + 0.25, 1.0 / ndim); 
            mds = 1;                                     
            if ((2 * ng - NDMX) >= 0)
            {
                mds = -1; 
                samples_2 = ng / NDMX + 1;
                nd = ng / samples_2;
                ng = samples_2 * nd;
            }
        }
        for (k = 1, i = 0; i < ndim; i++)
            k *= ng; 
        samples_2 = (samples/ k > 2) ? (samples/ k) : (2);
        calls = (double)samples_2 * (double)k;
        dxg = 1.0 / ng;
        for (dv2g = 1, i = 0; i < ndim; i++)
            dv2g *= dxg;
        dv2g = 1.0 / (samples_2 -1);
        xnd = nd;
        dxg *= xnd;
        xJac = 1.0 / calls;
        for (j = 0; j < ndim; j++)
        {
            dx[j] = regn[j + ndim] - regn[j];
            xJac *= dx[j];
        }
        
        if (nd != ndo)
        { 
            for (i = 0; i < (nd > ndo ? nd : ndo); i++)
                r[i] = 1.0;
            for (j = 0; j < ndim; j++)
                rebin(ndo / xnd, nd, r, xin, xi[j]);
            ndo = nd; 
            
        }
    }
    for (it = ittot; it <= itmx + ittot - 1; it++)
    {
        Ab.ti = Ab.tsi = 0.0;
        for (j = 0; j < ndim; j++)
        {
            kg[j] = 1;
            for (i = 0; i < nd; i++)
                d[i][j] = 0.0;
                
        }
        for (;;)
        {
            Ax.fb = 0.0; 
            Ax.f2b = 0.0;
            Ax.samples_2 = 0;
            for (k = 0; k < samples_2; k++)
            {
                wgt = xJac;
                for (j = 0; j < ndim; j++) 
                {
                    xrand = gfsr_rand(0, 1); 
                    xn = (kg[j] - xrand) * dxg + 1.0;   
                    index[j] = ((int)xn < NDMX) ? ((int)xn) : (NDMX);
                    index[j] = (index[j] > 1) ? (index[j]) : (1);
                    if (index[j] > 1) 
                    {
                        xo = xi[j][index[j] - 1] - xi[j][index[j] - 2]; 
                        rc = xi[j][index[j] - 2] + (xn - index[j]) * xo; 
                    }
                    else 
                    {
                        xo = xi[j][index[j] - 1]; 
                        rc = (xn - index[j]) * xo; 
                    }

                    x[j] = regn[j] + rc * dx[j];
                    wgt *= xo * xnd; 
                }
                fxn(x, &f); 
                f *= wgt; 
                Ax.f2 = f * f;
                Ax.fb += f; 
                Ax.f2b += Ax.f2;
                for (j = 0; j < ndim; j++)
                {
                    if (mds >= 0)
                    { 
                        d[index[j] - 1][j] += Ax.f2;
                    } 
                }
            } 
            Ax.f2b = sqrt(Ax.f2b * samples_2);
            Ax.f2b = (Ax.f2b - Ax.fb) * (Ax.f2b + Ax.fb); 
            if (Ax.f2b <= 0.0)
                Ax.f2b = TINY;
            Ab.ti += Ax.fb;   
            Ab.tsi += Ax.f2b; 
            if (mds < 0)
            {
                for (j = 0; j < ndim; j++)
                    d[index[j] - 1][j] += Ax.f2b; 
            }
            for (k = ndim - 1; k >= 0; k--)
            {
                kg[k] %= ng; 
                if (++kg[k] != 1)
                    break;
            }
            if (k < 0)
                break;
        } 
        Ab.tsi *= dv2g;
        Ai.Wgt = 1.0 / Ab.tsi;            
        Ai.sInt += Ai.Wgt * Ab.ti;         
        Ai.sWgt += Ai.Wgt;                 
        tgral = Ai.sInt / Ai.sWgt;         
        sd = sqrt(1.0 / Ai.sWgt);
        *sdin = sd;
        *tgralin = tgral;
        for (j = 0; j < ndim; j++)
        {
            xo = d[0][j];             
            xn = d[1][j];             
            d[0][j] = (xo + xn) / 2.0; 
            dt[j] = d[0][j];          
            for (i = 1; i < nd - 1; i++)
            {
                rc = xo + xn;
                xo = xn;
                xn = d[i + 1][j];
                d[i][j] = (rc + xn) / 3.0; 
                dt[j] += d[i][j];
            }
            d[nd - 1][j] = (xo + xn) / 2.0; 
            dt[j] += d[nd - 1][j];         
        }
        for (j = 0; j < ndim; j++)
        {
            rc = 0.0;
            for (i = 0; i < nd; i++)
            {
                dem1++;
                if (d[i][j] < TINY)
                    d[i][j] = TINY;
                r[i] = pow((1.0 - d[i][j] / dt[j]) / (log(dt[j]) - log(d[i][j])), ALPH); 
                rc += r[i];
            }
            rebin(rc / xnd, nd, r, xin, xi[j]); 
        }
    }
    
    ittot += itmx;
    
}
void testfun(double x[DIMENSION], double *f)
{
    *f = 1 / ((x[0] - 0.75) * (x[0] - 0.75) + 1e-6) + 1 / ((x[0] - 0.5) * (x[0] - 0.5) + 1e-6) + 1 / ((x[0] - 0.25) * (x[0] - 0.25) + 1e-6);
    // f[0] = 1/(x[0]*x[1]*x[2]*x[3]);
    //f[0] = 1/(x[0]);
    //f[0] = 1 / x[0] + 20 * (exp(-pow(10, 4) * (x[0] - 1) * (x[0] - 1))); // have a high peak
    // f[0] = sin(x[0])*cos(x[1]);
    // f[0]=  1 / ((x[0] - 0.25) * (x[0] - 0.25) + 1e-6) ;
    return;
}
int main()
{

    int i;
    double *estim;   /* estimators for integrals                     */
    double *std_dev; /* standard deviations                          */
    double reg[2 * DIMENSION]; /* integration domain                 */

    for (i = 0; i < DIMENSION; i++)
    {
        // reg[i] = 0.0000000001;
        // reg[i + DIMENSION] = 2.0;
        reg[i] = -4;
        reg[i + DIMENSION] =4.0;
    }

    vegas(reg, DIMENSION, testfun, 0, 10000, 10, estim, std_dev);
    // printf("Result: %g +/- %g\n", *estim, *std_dev);
    // // //randomNumberTest();
    // printf("\n============================\n");
    vegas(reg, DIMENSION, testfun, 0, 10000, 200, estim, std_dev);
    printf("Result: %g +/- %g\n", *estim, *std_dev);
    return 0;
}

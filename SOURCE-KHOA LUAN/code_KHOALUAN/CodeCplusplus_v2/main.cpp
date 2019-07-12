
#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include <stdio.h>
#include "vegas.h"
using namespace std;
double pi = 3.14159265359;

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
    // f là ước tính giá trị của hàm tính tích phân tại x 
    public:
    double f2;                 /* f^2 -> f bình phương                                 */
    double fb;                 /* Tổng f bên trong 1 bin. Mỗi bin được gieo 2 điểm     */
    double f2b;                /* Tổng f2 bên trong 1 bin.                             */
    unsigned long samples_2;   /* số điểm gieo bên trong 1 bin                         */
    
};
double DieuChinhDiem(double a,double b,double c,double d, double x){

    x = ((x-a)/(b-a))*(d-c) +c;
    return x;
}

void rebin(double wgt_avg, int nd, double wgt[], double xin[], double xi[])
{
    int i;
    int k = 0;
    double residues = 0.0;
    double xn = 0.0;
    double xo = 0.0;
    const char *filePath1 = "C:/Users/Ha Vy/Desktop/xi.txt";
    // const char *filePath2 = "C:/Users/Ha Vy/Desktop/success_test/xin.txt";
    FILE *file1;
    // FILE *file2;
    file1 = fopen(filePath1, "w+");
    // file2 = fopen(filePath2, "w+");
    for (i = 0; i < nd - 1; i++)
    {
        while (wgt_avg > residues)
        {
            residues += wgt[k++];
            // printf("\n========%lfw============\n", residues);
        }
        if (k > 1)
            xo = xi[k - 2];
        xn = xi[k - 1];
        residues -= wgt_avg; //residues of accummulate weights up to bin k-1
        xin[i] = xn - (xn - xo) * residues / wgt[k - 1];
        // printf("\n========%lf============\n", xin[i]);
    }
    for (i = 0; i < nd - 1; i++)
    {
        xi[i] = xin[i];
        // if(i == 0 || i == 20 || i == 40 || i == 60 || i == 80 )
        fprintf(file1, "%f\n", xin[i]);
    }
    xi[nd - 1] = 1.0;
    fprintf(file1, "%f\n", xin[nd - 1]);
    fclose(file1);
    // fclose(file2);
}
//=======================
void rebin1(double wgt_avg, int nd, double wgt[], double xin[], double xi[])
{
    int i;
    int k = 0;
    double residues = 0.0;
    double xn = 0.0;
    double xo = 0.0;
    const char *filePath1 = "C:/Users/Ha Vy/Desktop/xi_1.txt";
    // const char *filePath2 = "C:/Users/Ha Vy/Desktop/success_test/xin.txt";
    FILE *file1;
    // FILE *file2;
    file1 = fopen(filePath1, "w+");
    // file2 = fopen(filePath2, "w+");
    for (i = 0; i < nd - 1; i++)
    {
        while (wgt_avg > residues)
        {
            residues += wgt[k++];
            // printf("\n========%lfw============\n", residues);
        }
        if (k > 1)
            xo = xi[k - 2];
        xn = xi[k - 1];
        residues -= wgt_avg; //residues of accummulate weights up to bin k-1
        xin[i] = xn - (xn - xo) * residues / wgt[k - 1];
        // printf("\n========%lf============\n", xin[i]);
    }
    for (i = 0; i < nd - 1; i++)
    {
        xi[i] = xin[i];
        // if(i == 0 || i == 20 || i == 40 || i == 60 || i == 80 )
        fprintf(file1, "%f\n", xin[i]);
    }
    xi[nd - 1] = 1.0;
    fprintf(file1, "%f\n", xin[nd - 1]);
    fclose(file1);
    // fclose(file2);
}
//=============================
void rebin2(double wgt_avg, int nd, double wgt[], double xin[], double xi[])
{
    int i;
    int k = 0;
    double residues = 0.0;
    double xn = 0.0;
    double xo = 0.0;
    const char *filePath1 = "C:/Users/Ha Vy/Desktop/xi_2.txt";
    // const char *filePath2 = "C:/Users/Ha Vy/Desktop/success_test/xin.txt";
    FILE *file1;
    // FILE *file2;
    file1 = fopen(filePath1, "w+");
    // file2 = fopen(filePath2, "w+");
    for (i = 0; i < nd - 1; i++)
    {
        while (wgt_avg > residues)
        {
            residues += wgt[k++];
            // printf("\n========%lfw============\n", residues);
        }
        if (k > 1)
            xo = xi[k - 2];
        xn = xi[k - 1];
        residues -= wgt_avg; //residues of accummulate weights up to bin k-1
        xin[i] = xn - (xn - xo) * residues / wgt[k - 1];
        // printf("\n========%lf============\n", xin[i]);
    }
    for (i = 0; i < nd - 1; i++)
    {
        xi[i] = xin[i];
        // if(i == 0 || i == 20 || i == 40 || i == 60 || i == 80 )
        fprintf(file1, "%f\n", xin[i]);
    }
    xi[nd - 1] = 1.0;
    fprintf(file1, "%f\n", xin[nd - 1]);
    fclose(file1);
    // fclose(file2);
}
//=============================
void rebin3(double wgt_avg, int nd, double wgt[], double xin[], double xi[])
{
    int i;
    int k = 0;
    double residues = 0.0;
    double xn = 0.0;
    double xo = 0.0;
    const char *filePath1 = "C:/Users/Ha Vy/Desktop/xi_3.txt";
    // const char *filePath2 = "C:/Users/Ha Vy/Desktop/success_test/xin.txt";
    FILE *file1;
    // FILE *file2;
    file1 = fopen(filePath1, "w+");
    // file2 = fopen(filePath2, "w+");
    for (i = 0; i < nd - 1; i++)
    {
        while (wgt_avg > residues)
        {
            residues += wgt[k++];
            // printf("\n========%lfw============\n", residues);
        }
        if (k > 1)
            xo = xi[k - 2];
        xn = xi[k - 1];
        residues -= wgt_avg; //residues of accummulate weights up to bin k-1
        xin[i] = xn - (xn - xo) * residues / wgt[k - 1];
        // printf("\n========%lf============\n", xin[i]);
    }
    for (i = 0; i < nd - 1; i++)
    {
        xi[i] = xin[i];
        // if(i == 0 || i == 20 || i == 40 || i == 60 || i == 80 )
        fprintf(file1, "%f\n", xin[i]);
    }
    xi[nd - 1] = 1.0;
    fprintf(file1, "%f\n", xin[nd - 1]);
    fclose(file1);
    // fclose(file2);
}
 double gfsr_rand (double a, double b){

   return a + (b - a)*rand()/RAND_MAX;

 }
void vegas(double regn[], int ndim, void (*fxn)(double x[], double *f), //replace f[] by f(1)
           int init, unsigned long samples, int itmx,          // fxn(x, &f); /* call integrand at point x */
           double *tgralin, double *sdin)
{

    double tgral, sd;
    static int ndo;   /*                                    (ndo)     */
    int it;           /* iteration counter                  (it)      */
    static int ittot; /* iteration counter across init>1              */
    int i, j, k;      /* counters                           (i, j, k) */
    int nd;           /* slices in grid (c.f. NDMX)         (nd)      */
    int ng;           /* real number of bins on an axis     (ng)      */
    unsigned int samples_2; /* number of calls within bin         (samples_2)     */
    static int mds;   /* ==1: statified smpl.               (mds)     */
    int index[MXDIM];    /* index of bin (ng bin per axis)     (index[])    */
    int kg[MXDIM];    /* contral samples_2 calls locate same bin  (kg[])    */
    double calls;     /* real total number of calls to fxn  (calls)   */
    double dv2g;      /*                                    (dv2g)    */
    double dxg;       /* ratio of number of grids per bins  (dxg)     */
    double rc,rc1;        /*                                    (rc)      */
    double wgt;       /* weight                             (wgt)     */
    double xn,xn1;        /*                                    (xn)      */
    double xnd;       /* another name of nd                 (xnd)     */
    double xo,xo1;        /*                                    (xo)      */
    double xJac;      /* Jacobian of integration            (xjac)    */
    static info_iteration Ai; /* ...for integrand                             */
    static info_bins Ab;     /* ...for integrand                             */
    static info_bin Ax;                  /* ...for integrand                             */
    double f;                      /* passed into fxn for evaluation at x          */
    double x[MXDIM];               /* evaluation point                   (x[])     */
    double d[NDMX][MXDIM];         /*                                    (d[][])   */
    // double di[NDMX][MXDIM];        /* delta i                            (di[][])  */
    double dt[MXDIM],dit[MXDIM];              /*                                    (dt[])    */
    double r[NDMX];                /*                                    (r[])     */
    static double xi[MXDIM][NDMX]; /*                                    (xi[][])  */
    double xin[NDMX];              /* aux. variable for rebinning        (xin[])   */
    double dx[MXDIM];              /* width of integration region        (dx[])    */
    double xrand;                  /* uniform dens 0.0 <= xrand < 1.0   (xrand)   */

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
        ittot = 1;
        // printf("taolao");
    }
    if (init <= 2)
    { /* inherit grid and results    */
        nd = NDMX;
        ng = 1;
        if (mds)
        {
            ng = (int)pow(samples/ 2.0 + 0.25, 1.0 / ndim); // M=2N^n; ng <-- N
            mds = 1;                                       //concentrate subintervals where integrand is largest in magnitude
            if ((2 * ng - NDMX) >= 0)
            {
                mds = -1; //concentrate subintervals where integrand is largest in error (sigma)
                samples_2 = ng / NDMX + 1;
                nd = ng / samples_2;
                ng = samples_2 * nd;
            }
        }
        for (k = 1, i = 0; i < ndim; i++)
            k *= ng; // k <-- ng^ndim
        samples_2 = (samples/ k > 2) ? (samples/ k) : (2);
        calls = (double)samples_2 * (double)k;
        dxg = 1.0 / ng;
        for (dv2g = 1, i = 0; i < ndim; i++)
            dv2g *= dxg;
        dv2g = 1.0 / (samples_2 -1);
        xnd = nd;
        dxg *= xnd; // dxg <-- nd/ng
        xJac = 1.0 / calls;
        for (j = 0; j < ndim; j++)
        {
            dx[j] = regn[j + ndim] - regn[j];
            xJac *= dx[j];
        }
        // printf("\nndo _1 = %d\n", ndo);
        // printf("\nndo _1 = %lf\n", xnd);
        if (nd != ndo)
        { //ndo=1 at this step: initialize bin
            for (i = 0; i < (nd > ndo ? nd : ndo); i++)
                r[i] = 1.0;
            for (j = 0; j < ndim; j++)
                rebin(ndo / xnd, nd, r, xin, xi[j]);
            ndo = nd; //number of subdivisions on an axis
            // printf("\nndo _1 = %d\n", ndo);
            // printf("\nnxd _1 = %lf\n", xnd);
            // printf("\nnd _1 = %d\n", nd);
        }
    }
    // printf("\nnd = %d\n", nd);
    // printf("\nng = %d\n", ng);
    // printf("\nxnd = %lf\n", xnd);
     const char *filePath3 = "C:/Users/Ha Vy/Desktop/index.txt";
            FILE *file3;
            file3 = fopen(filePath3, "w+");
    for (it = ittot; it <= itmx + ittot - 1; it++)
    {
        Ab.ti = Ab.tsi = 0.0;
        for (j = 0; j < ndim; j++)
        {
            kg[j] = 1;
            for (i = 0; i < nd; i++)
                d[i][j] = 0.0;
                // di[i][j] = 0.0;
        }
        for (;;)
        {
            Ax.fb = 0.0; //refresh to start a hypercube
            Ax.f2b = 0.0;
            Ax.samples_2 = 0;
            // printf("\nkg[j] = %d\n",kg[0]);
            for (k = 0; k < samples_2; k++)
            {
                wgt = xJac; // tạo điểm ngẫu nhiên X 
                for (j = 0; j < ndim; j++) 
                {
                    xrand = gfsr_rand(0, 1); //uniform random number, prob <--1/ng
                    xn = (kg[j] - xrand) * dxg + 1.0;   //1<kg[j]<=ng // xn: số chay theo ng (bởi kg theo ng)
                    // printf("\n%d\n",kg[j]);
                    // // // xn: là cách bâm đoạn nhỏ nên 1<=xn<NDMX
                    // // // index[dim] được gán bằng xn ==> idsu là chỉ số để chạy trên trục N 
                    if(k==0 && j==0) fprintf(file3, "%d\n", index[j] );;
                    index[j] = ((int)xn < NDMX) ? ((int)xn) : (NDMX);
                    index[j] = (index[j] > 1) ? (index[j]) : (1);
                    if (index[j] > 1) // chỉ số khac 1 (đang nằm trong [a,b])
                    {
                        xo = xi[j][index[j] - 1] - xi[j][index[j] - 2];  // xo = xi[i+1] - xi[i]
                        rc = xi[j][index[j] - 2] + (xn - index[j]) * xo; //
                    }
                    else // bắt đâu [a,b]
                    {
                        xo = xi[j][index[j] - 1]; // => xo = xi[dim][0]
                        rc = (xn - index[j]) * xo; 
                    }
                    // x[j] = DieuChinhDiem(0.0, 1.0, xi[j][index[j] - 1], xi[j][index[j] - 2],xrand);
                    x[j] = regn[j] + rc * dx[j];
                    wgt *= xo * xnd; // p(x)=1/Ndelt_i trong cong thức ??? p(xrand)dxrand = p(rc)drc  ???, where p(xrand) = 1/ng
                    //wgt nhân dồn cho các chiều tương tự như 1 chiều
                }
                // tạo điểm ngẫu nhiên X
                fxn(x, &f); /* call integrand at point x */
                // if (f != 0.0) // nếu f ==0 thì gieo điểm lại
                //     ++Ax.samples_2;
                f *= wgt; // contribution of x to I // f ở đây tương đương với việc giá trị f nhân với 1 xác suất
                Ax.f2 = f * f;
                Ax.fb += f; //sum over samples_2 points whithin hypercube-->contribution of an hypercube to I
                Ax.f2b += Ax.f2;
                // di[index_Ô][DIM] : chứa Sum f
                // d[index_Ô][DIM] : chứa Sum f2
                for (j = 0; j < ndim; j++)
                {
                    // di[index[j] - 1][j] += f; // cộng dồn là do 2 điểm gieo
                    if (mds >= 0)
                    { //printf("\nKhong c\n");
                        d[index[j] - 1][j] += Ax.f2;
                    } //sum up samples_2 magnitude f2 whithin each bin of a hypercube
                }
            } /* end of loop within hypercube */
            // printf("\nGieo diem thu %d\n", k);
            // Ax.f2b = (Ax.fb - Ax.f2b)/(samples_2-1);
            Ax.f2b = sqrt(Ax.f2b * samples_2);
            Ax.f2b = (Ax.f2b - Ax.fb) * (Ax.f2b + Ax.fb); //sum{aCube}{f^2}-(sum{aCube}f)^2-->  acube
            if (Ax.f2b <= 0.0)
                Ax.f2b = TINY;
            Ab.ti += Ax.fb;   //sum over hypercubes --> sum{over calls points}f --> I of an iteration
            Ab.tsi += Ax.f2b; //sum(cubes){sigma^2 of cube) --> sigtot^2 of an iteration
            if (mds < 0)
            {
                for (j = 0; j < ndim; j++)
                    d[index[j] - 1][j] += Ax.f2b; //sum(over cubs per bin){sigma^2)-->sigbin^2
            }
            for (k = ndim - 1; k >= 0; k--)
            {
                kg[k] %= ng; //create kg[k]: 1<=kg[k]<=ng; ng^ndim values (equal or none)
                // printf("\n%d\n",kg[k]);
                if (++kg[k] != 1)
                    break;
            }
            if (k < 0)
                break;
        } /* end of loop over hypercubes */
        fclose(file3);
        // printf("\nng == %d\n",ng);
        int dem1=0;
        Ab.tsi *= dv2g;
        Ai.Wgt = 1.0 / Ab.tsi;             // ~ 1/(sigtot_i)^2 of iteration ith
        Ai.sInt += Ai.Wgt * Ab.ti;         // ~ sum(iteration){I_i/(sigtot_i)^2}
        Ai.sWgt += Ai.Wgt;                 // ~ sum(iteration){1/(sigtot_i)^2}-->1/sigI^2
        tgral = Ai.sInt / Ai.sWgt;         //~sigI^2.sum(iteration){I_i/(sigtot_i)^2}
        sd = sqrt(1.0 / Ai.sWgt); //sigI^2 // độ lệch chuẩn 
        *sdin = sd;
        *tgralin = tgral;
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
                dem1++;
                if (d[i][j] < TINY)
                    d[i][j] = TINY;
                // r[i] = d[i][j] / dt[j];
                r[i] = pow((1.0 - d[i][j] / dt[j]) / (log(dt[j]) - log(d[i][j])), ALPH); //sigbi contribution of bin ith in sumsigbin
                rc += r[i];
                
            }
            // printf("\nrc = %lf\n", rc);
            rebin(rc / xnd, nd, r, xin, xi[j]); //rebin wrt input ratio sigbi/
            
            switch (j)
            {
            case /* constant-expression */0:
                /* code */rebin(rc / xnd, nd, r, xin, xi[j]);
                break;
            case /* constant-expression */1:
                /* code */rebin1(rc / xnd, nd, r, xin, xi[j]);
                break;
            case /* constant-expression */2:
                /* code */rebin2(rc / xnd, nd, r, xin, xi[j]);
                break;
            case /* constant-expression */3:
                /* code */rebin3(rc / xnd, nd, r, xin, xi[j]);
                break;
            default:
                break;
            }
        }
    }
    
    ittot += itmx;
    
}
void testfun(double x[DIMENSION], double f[FUNCTIONS])
{
    // f[0] = 1 / ((x[0] - 0.75) * (x[0] - 0.75) + 1e-6) + 1 / ((x[0] - 0.5) * (x[0] - 0.5) + 1e-6) + 1 / ((x[0] - 0.25) * (x[0] - 0.25) + 1e-6);
    f[0] = 1/(x[0]*x[1]*x[2]+1e-6);
    //f[0] = 1/(x[0]);
    // f[0] = 1 / x[0] + 20 * (exp(-pow(10, 4) * (x[0] - 1) * (x[0] - 1))); // have a high peak
    // f[0] = sin(x[0])*cos(x[1]);
    // f[0]=  1 / ((x[0] - 0.25) * (x[0] - 0.25) + 1e-6) ;
    //2d
    f[0]=exp(-20*(x[0]*x[0]+x[1]*x[1]));
    //=======4d
    // double a;
    // a = pow(1/(0.1*pow(pi,0.5)),4);
    // f[0] = a*exp(-pow(0.1,-2)*((x[0]-0.5)*(x[0]-0.5)+(x[1]-0.5)*(x[1]-0.5)+(x[2]-0.5)*(x[2]-0.5)+(x[3]-0.5)*(x[3]-0.5))) ;
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
        reg[i] = -1.0;
        reg[i + DIMENSION] = 1.0;
        // reg[i] = 0.0;
        // reg[i + DIMENSION] =1.0;
    }

//==========4dim gaussian
   
    vegas(reg, DIMENSION, testfun, 0, 10000, 10, estim, std_dev);
    // printf("Result: %g +/- %g\n", *estim, *std_dev);
    // // //randomNumberTest();
    // printf("\n============================\n");
    vegas(reg, DIMENSION, testfun, 3, 10000, 90, estim, std_dev);
    printf("Result: %g +/- %g\n", *estim, *std_dev);
    return 0;
}

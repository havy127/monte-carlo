Vegas.h
NDMX 100 : số đoạn trên 1 trục

================================================
Tham số của hàm Vegas

parameter[1]: regn[] 
parameter[2]: ndim
parameter[3]: void (*fxn)(double x[], double *f)
parameter[4]: init
parameter[5]: ncall                            == 10^4             | Tổng số mẫu M
parameter[6]: itmx
parameter[7]:  nprn ----------- NPRN_INPUT | NPRN_RESULT == 3
parameter[8]: *tgralin
parameter[9]: *sdin
parameter[10]: *chi2ain

============================================
Biến bên trong hàm Vegas

static int ittot: vòng lặp
static int mds: ==1: statified smpl.
static int ndo:

double tgral, sd, chi2a
*tgralin = tgral 
*sdin = sd
*chi2ain = chi2a
ng: tổng số Ô được chia, mỗi bin gieo 2 điểm : M=2N^n; ng <-- N
npg: số điểm gieo vào 1 Ô
ia[MXDIM]:  index of bin (ng bin per axis)
kg[MXDIM]: contral npg calls locate same bin
calls: Tổng số điểm M sau khi tính toán . Bị hao hụt
    tổng số calls to fxn (tổng số điểm gieo : bins * 2)
dv2g: ????????
dxg: 1/ng (1/ số đoạn lớn của mỗi chiều)
    tương đương trong công thức 1/N 
rc: chứa Tổng mi của từng trục
wgt: 
xn: 
nd: ???????????????????
    nd = NDMX 
xnd: another name of nd
xo
xJac: Jacobian of integration 
f: kết quả integrand tại x 
x[]
d[][]
di[][]
dt[]: mãng chưa fitb
r[]: mãng chứa các mi của từng trục
xi[][]
xin[]
dx[]

*******STRUCT*********
1. iterAccu Ai {

    wgt
    sWgt: cumulative sum for weights
    sChi: cumulative sum for chi^2
    sInt: cumulative sum for integral

} // thông số tích lũy ở cuối vòng lặp

2. binAccu Ab{
    ti: tổng f qua Tổng các bins
    tsi: tổng variances các bins

} thông số tông qua các bins 

3. pointAccu Ax{

    f2: f bình phương
    fb: tổng f bên trong 1 bin 
    f2b: tổng f^2 bên trong 1 bin 
    npg: số mẫu trong 1 bin Nếu f !=0

} thông số bên trong các bins
=============================================
Tham số của hàm rebin 

[1]:wgt_avg
[2]:nd  
[3]:wgt[]
[4]:xin[]
[5]:xi[]

Biến bên trong hàm rebin
i : biến đếm vòng for
k: 
xo: 
xn:
residues: ???? 


===================================
Dòng 181
k: tổng số Ô của lưới được tính theo công thức k <--- ng^ndim
Dòng 196
xJac: là trọng số của 1 hạt ==>xJac = 1/(calls*(b-a)*(d-c)*....)
Dòng 357:
kg[ndim]: chạy theo ndim -> kg[0]:kg[ndim] <- 1

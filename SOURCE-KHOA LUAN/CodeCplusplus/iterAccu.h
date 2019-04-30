class iterAccu         /* accumulator for stuff at end of iteration... */
{
private:
    /* data */
public:
    double Wgt;       /* weight                             (wgt)     */
    double sWgt;      /* cumulative sum for weights         (swgt)    */
    double sChi;      /* cumulative sum for chi^2           (schi)    */
    double sInt;      /* cumulative sum for integral        (si)      */
};


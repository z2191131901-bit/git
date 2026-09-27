/* 普通最小二乘 y=k*x+b：假定x准确、误差主要在y。
 * 在线累计中心化二阶量，避免 sum(x*x)-sum(x)^2/n 的大数相减。 */
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
static int fit(const double *x,const double *y,size_t n,double *slope,double *offset)
{
    if(!x || !y || !slope || !offset || n<2) return 0;
    double mx=0,my=0,sxx=0,sxy=0;
    for(size_t i=0;i<n;i++) {
        if(!isfinite(x[i]) || !isfinite(y[i])) return 0;
        double dx=x[i]-mx,dy=y[i]-my;
        mx+=dx/(double)(i+1); my+=dy/(double)(i+1);
        sxx+=dx*(x[i]-mx);
        sxy+=dx*(y[i]-my);
    }
    /* 所有x相同，斜率无法确定；极小尺度也拒绝，避免除以近零。 */
    if(!isfinite(sxx) || !isfinite(sxy) || sxx<=DBL_MIN) return 0;
    double k=sxy/sxx,b=my-k*mx;
    if(!isfinite(k) || !isfinite(b)) return 0;
    *slope=k; *offset=b;
    return 1;
}
int main(void)
{
    const double x[]={0,1,2,3},y[]={1,3,5,7},same[]={1,1,1,1};
    double k=0,b=0;
    assert(fit(x,y,4,&k,&b) && k==2 && b==1);
    const double noisy[]={1.1,2.9,4.9,7.1};
    assert(fit(x,noisy,4,&k,&b) && fabs(k-2)<1e-12 && fabs(b-1)<1e-12);
    assert(!fit(same,y,4,&k,&b) && fabs(k-2)<1e-12);
    const double shifted[]={1e9,1e9+1,1e9+2,1e9+3};
    assert(fit(shifted,y,4,&k,&b) && k==2 && b==1-2e9);
    assert(!fit(x,y,1,&k,&b));
    puts("least squares: PASS, slope=2 intercept=1");
    return 0;
}

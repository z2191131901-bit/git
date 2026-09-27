/* 双积分对象的滑模控制：s=v+lambda*(p-r)，固定位置参考。
 * 比较没有切换项、符号切换、饱和边界层三种形式。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
typedef struct { double error,variation; } Result;
static int control(double position,double velocity,double reference,
                   double eta,double width,double *output)
{
    if(!output || !isfinite(position) || !isfinite(velocity) ||
       !isfinite(reference) || !isfinite(eta) || eta<0 ||
       !isfinite(width) || width<0) return 0;
    const double lambda=2,k=3;
    double s=velocity+lambda*(position-reference);
    /* width=0用sign；非零边界层内用s/width，层外饱和为±1。 */
    double switching=width==0 ? (s>0 ? 1 : (s<0 ? -1 : 0))
                             : fmax(-1,fmin(1,s/width));
    double u=-lambda*velocity-k*s-eta*switching;
    if(!isfinite(u)) return 0;
    *output=u; return 1;
}
static Result run(double eta,double width)
{
    double p=0,v=0,previous_u=0,total_variation=0;
    const double dt=.001;
    for(int k=0;k<15000;k++) {
        double u;
        assert(control(p,v,1,eta,width,&u) && fabs(u)<20);
        double acceleration=u+(k>=2000 ? .3 : 0);
        p+=v*dt+acceleration*dt*dt/2; v+=acceleration*dt;
        /* 只统计最后5秒，排除正常启动过程。没有测量噪声。 */
        if(k>=10000) total_variation+=fabs(u-previous_u);
        previous_u=u;
    }
    Result result={fabs(p-1),total_variation}; return result;
}
int main(void)
{
    double u=99;
    assert(control(1.05,0,1,.8,.2,&u));
    assert(fabs(u+.7)<1e-12); /* s=.1，-3*.1-.8*(.1/.2)=-.7。 */
    assert(control(1,0,1,.8,0,&u) && u==0);
    assert(!control(0,0,1,.8,-1,&u) && u==0);
    Result nominal=run(0,0),sign=run(.8,0),smooth=run(.8,.05);
    assert(fabs(nominal.error-.05)<1e-5);
    assert(smooth.error<nominal.error*.3 && sign.error<.002);
    assert(fabs(smooth.error-.3/(2*(3+.8/.05)))<1e-5);
    assert(sign.variation>10*smooth.variation);
    printf("SMC: final error nominal=%.6f sign=%.6f boundary=%.6f\n",
           nominal.error,sign.error,smooth.error);
    printf("SMC: last 5s control variation sign=%.6f boundary=%.6f\n",
           sign.variation,smooth.variation);
    return 0;
}

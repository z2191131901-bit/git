/* 大误差时暂停误差积分，小误差时恢复；暂停不等于清零。
 * 为公平比较，两组控制都使用相同的条件积分抗饱和和测量微分。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
typedef struct { double integral,previous,derivative; int initialized; } Pid;
static double clamp(double x,double lo,double hi) { return fmax(lo,fmin(hi,x)); }
static int control(Pid *s,double reference,double y,double dt,
                   double threshold,int separation,double *output)
{
    const double kp=2,ki=4,kd=.1,tau=.03,limit=2;
    if(!s || !output || !isfinite(reference) || !isfinite(y) ||
       !isfinite(dt) || dt<=0 || !isfinite(threshold) || threshold<0) return 0;
    Pid next=*s;
    double e=reference-y;
    double rate=s->initialized ? (y-s->previous)/dt : 0;
    next.derivative=s->initialized
        ? s->derivative+dt/(tau+dt)*(rate-s->derivative) : 0;
    /* 门限单位与误差相同。误差超限仅冻结积分贡献，保留已有负载补偿。 */
    double increment=(!separation || fabs(e)<=threshold) ? ki*e*dt : 0;
    next.integral=s->integral+increment;
    double raw=kp*e+next.integral-kd*next.derivative;
    if((raw>limit && increment>0) || (raw<-limit && increment<0)) {
        next.integral=s->integral;
        raw=kp*e+next.integral-kd*next.derivative;
    }
    if(!isfinite(raw) || !isfinite(next.integral)) return 0;
    next.previous=y; next.initialized=1;
    *s=next; *output=clamp(raw,-limit,limit);
    return 1;
}
static double run(int separate,double *final)
{
    Pid pid={0};
    double y=0,peak=0;
    const double dt=.005,a=exp(-dt);
    for(int k=0;k<3000;k++) {
        double u;
        assert(control(&pid,1,y,dt,.4,separate,&u));
        assert(fabs(u)<=2);
        y=a*y+(1-a)*u; /* y'=-y+u，区间内输入恒定的精确离散式。 */
        peak=fmax(peak,y);
    }
    *final=y;
    return peak;
}
int main(void)
{
    Pid p={.2,0,0,0};
    double u=99;
    assert(control(&p,1,0,.01,.4,1,&u) && p.integral==.2);
    assert(control(&p,.1,0,.01,.4,1,&u) && fabs(p.integral-.204)<1e-12);
    double before=p.integral;
    assert(!control(&p,1,0,0,.4,1,&u) && p.integral==before);
    double normal_end,separate_end;
    double normal_peak=run(0,&normal_end),separate_peak=run(1,&separate_end);
    assert(fabs(normal_end-1)<.001 && fabs(separate_end-1)<.001);
    assert(separate_peak<normal_peak);
    printf("integral separation: normal peak=%.6f separated peak=%.6f final=%.6f\n",
           normal_peak,separate_peak,separate_end);
    return 0;
}

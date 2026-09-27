/* 变化率限制：每秒最大上升/下降速度，可使用不等时间间隔。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
static int slew(double *value,double target,double rise,double fall,double dt)
{
    if(!value || !isfinite(*value) || !isfinite(target) ||
       !isfinite(rise) || !isfinite(fall) || !isfinite(dt) ||
       rise<0 || fall<0 || dt<=0) return 0;
    double delta=target-*value, up=rise*dt, down=fall*dt;
    if(!isfinite(delta) || !isfinite(up) || !isfinite(down)) return 0;
    /* 只限制本次允许的变化量；最后一步自动缩短，不越过目标。 */
    delta=fmax(-down,fmin(up,delta));
    double next=*value+delta;
    if(!isfinite(next)) return 0;
    *value=next;
    return 1;
}
int main(void)
{
    double y=0;
    assert(slew(&y,1,2,4,.1) && fabs(y-.2)<1e-12);
    assert(slew(&y,1,2,4,.4) && fabs(y-1)<1e-12);
    assert(slew(&y,0,2,4,.1) && fabs(y-.6)<1e-12);
    assert(slew(&y,.5,2,4,.1) && y==.5);
    assert(!slew(&y,0,2,4,0) && y==.5);
    assert(slew(&y,1,0,4,.1) && y==.5);
    puts("slew: PASS, asymmetric rates and no overshoot");
    return 0;
}

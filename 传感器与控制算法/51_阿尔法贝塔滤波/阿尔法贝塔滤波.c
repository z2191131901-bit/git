/* 常速度模型的固定增益位置/速度估计；不是自适应卡尔曼滤波。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
typedef struct { double position,velocity; } Estimate;
static int step(Estimate *s,double dt,int has_measurement,double measured,
                double alpha,double beta)
{
    if(!s || !isfinite(s->position) || !isfinite(s->velocity) ||
       !isfinite(dt) || dt<=0 || !isfinite(alpha) || !isfinite(beta) ||
       alpha<=0 || alpha>=1 || beta<=0 || beta>=4-2*alpha ||
       (has_measurement!=0 && has_measurement!=1) ||
       (has_measurement && !isfinite(measured))) return 0;
    double predicted=s->position+s->velocity*dt;
    Estimate next={predicted,s->velocity};
    if(has_measurement) {
        double residual=measured-predicted;
        next.position+=alpha*residual;
        next.velocity+=beta*residual/dt;
    }
    if(!isfinite(next.position) || !isfinite(next.velocity)) return 0;
    *s=next;
    return 1;
}
int main(void)
{
    Estimate s={0,1};
    assert(step(&s,.1,1,.2,.5,.1));
    assert(fabs(s.position-.15)<1e-12 && fabs(s.velocity-1.1)<1e-12);
    assert(step(&s,.1,0,NAN,.5,.1)); /* 缺测只预测，不把0当成位置。 */
    assert(fabs(s.position-.26)<1e-12);
    assert(!step(&s,0,1,0,.5,.1));
    Estimate track={0,0};
    for(int i=1;i<=200;i++) assert(step(&track,.1,1,2*i*.1,.5,.1));
    assert(fabs(track.position-40)<1e-8 && fabs(track.velocity-2)<1e-8);
    puts("alpha-beta: PASS, constant-speed estimate=2 m/s");
    return 0;
}

/* 手动时用实际手动输出反推积分；切入自动当拍对齐上次施加量。
 * 自动后继续正常PID。无扰指切换当拍输出连续，不代表后续输出不变化。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
typedef struct {
    double integral,previous,derivative,last_output;
    int initialized,was_auto;
} Controller;
static double clamp(double x) { return fmax(-1,fmin(1,x)); }
static int control(Controller *s,double r,double y,double manual,int automatic,
                   int bumpless,double dt,double *output)
{
    const double kp=2,ki=1,kd=.1,tau=.05;
    if(!s || !output || !isfinite(r) || !isfinite(y) || !isfinite(manual) ||
       !isfinite(dt) || dt<=0 || (automatic!=0 && automatic!=1)) return 0;
    Controller next=*s;
    double e=r-y,rate=s->initialized ? (y-s->previous)/dt : 0;
    next.derivative=s->initialized ? s->derivative+dt/(tau+dt)*(rate-s->derivative) : 0;
    double pd=kp*e-kd*next.derivative,u;
    if(!automatic) {
        u=clamp(manual);
        /* 将I对齐到实际限幅后的手动输出，而不是未限幅的手动命令。 */
        if(bumpless) next.integral=u-pd;
    } else if(bumpless && s->initialized && !s->was_auto) {
        /* 自动接管当拍不累加误差积分，严格保持上次实际输出。 */
        u=s->last_output;
        next.integral=u-pd;
    } else {
        double increment=ki*e*dt;
        next.integral+=increment;
        double raw=pd+next.integral;
        if((raw>1 && increment>0) || (raw<-1 && increment<0))
            next.integral=s->integral;
        u=clamp(pd+next.integral);
    }
    if(!isfinite(pd) || !isfinite(next.integral) || !isfinite(u)) return 0;
    next.previous=y; next.initialized=1; next.was_auto=automatic; next.last_output=u;
    *s=next; *output=u;
    return 1;
}
static double run(int bumpless,double *final)
{
    Controller s={0}; double y=0,jump=0;
    const double dt=.01,a=exp(-dt);
    for(int k=0;k<2500;k++) {
        double before=s.last_output,u;
        assert(control(&s,.7,y,.4,k>=200,bumpless,dt,&u));
        if(k==200) jump=fabs(u-before);
        y=a*y+(1-a)*u;
    }
    *final=y;
    return jump;
}
int main(void)
{
    Controller s={0}; double u;
    assert(control(&s,1,.4,.4,0,1,.01,&u) && u==.4);
    assert(fabs(s.integral+.8)<1e-12);
    assert(control(&s,1,.4,.4,1,1,.01,&u) && u==.4);
    assert(control(&s,1,.4,.4,1,1,.01,&u) && fabs(u-.406)<1e-12);
    assert(!control(&s,1,.4,.4,1,1,0,&u));
    double end1,end2,jump1=run(0,&end1),jump2=run(1,&end2);
    assert(jump1>.1 && jump2<1e-12);
    assert(fabs(end1-.7)<.001 && fabs(end2-.7)<.001);
    printf("bumpless transfer: direct jump=%.6f tracked jump=%.6f final=%.6f\n",
           jump1,jump2,end2);
    return 0;
}

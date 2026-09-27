/* I_next=I+dt*(Ki*e+Kb*(u_applied-u_raw))。
 * 本例实际执行器只有对称幅值限幅；若还有后级限速等，应反馈真实施加量。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
typedef struct { double integral,previous,derivative; int initialized; } Pid;
static int control(Pid *s,double reference,double y,double dt,double kb,double *output)
{
    const double kp=2,ki=1,kd=.05,tau=.05;
    if(!s || !output || !isfinite(reference) || !isfinite(y) ||
       !isfinite(dt) || dt<=0 || !isfinite(kb) || kb<0 || kb*dt>1) return 0;
    Pid next=*s;
    double e=reference-y;
    double rate=s->initialized ? (y-s->previous)/dt : 0;
    next.derivative=s->initialized ? s->derivative+dt/(tau+dt)*(rate-s->derivative) : 0;
    /* 先按旧积分算当前输出，再更新下一拍积分，避免隐含代数环。 */
    double raw=kp*e+s->integral-kd*next.derivative;
    double applied=fmax(-1,fmin(1,raw));
    next.integral=s->integral+dt*(ki*e+kb*(applied-raw));
    if(!isfinite(raw) || !isfinite(next.integral)) return 0;
    next.previous=y; next.initialized=1;
    *s=next; *output=applied;
    return 1;
}
static double run(double kb,double *integral_at_switch,double *final)
{
    Pid s={0};
    const double dt=.01,a=exp(-dt);
    double y=0,error_area=0;
    for(int k=0;k<2000;k++) {
        double r=k<500 ? 2 : .5,u;
        if(k==500) *integral_at_switch=s.integral;
        assert(control(&s,r,y,dt,kb,&u));
        assert(fabs(u)<=1);
        y=a*y+(1-a)*u;
        if(k>=500 && k<1000) error_area+=fabs(r-y)*dt;
    }
    *final=y;
    return error_area;
}
int main(void)
{
    Pid s={0};
    double u;
    /* raw=4, applied=1，I_next=.01*(2+2*(1-4))=-.04。 */
    assert(control(&s,2,0,.01,2,&u) && u==1 && fabs(s.integral+.04)<1e-12);
    assert(!control(&s,2,0,.01,200,&u));
    double old_i,new_i,old_end,new_end;
    double old_area=run(0,&old_i,&old_end),new_area=run(2,&new_i,&new_end);
    assert(fabs(new_i)<fabs(old_i) && new_area<old_area*.5);
    assert(fabs(new_end-.5)<.01);
    printf("back calculation: switch I without=%.6f with=%.6f; recovery IAE=%.6f / %.6f\n",
           old_i,new_i,old_area,new_area);
    return 0;
}

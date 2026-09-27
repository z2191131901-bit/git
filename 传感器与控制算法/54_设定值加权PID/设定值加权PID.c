/* 两自由度PID的常用形式：P=Kp*(b*r-y)，I仍积完整误差，D只对测量。
 * 微分设定值权重固定为0；比例设定值权重b可调。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
typedef struct { double integral,previous,derivative; int initialized; } Pid;
static int control(Pid *s,double r,double y,double dt,double weight,double *output)
{
    const double kp=4,ki=2,kd=1,tau=.02;
    if(!s || !output || !isfinite(r) || !isfinite(y) || !isfinite(dt) ||
       dt<=0 || !isfinite(weight) || weight<0 || weight>1) return 0;
    Pid next=*s;
    double rate=s->initialized ? (y-s->previous)/dt : 0;
    next.derivative=s->initialized ? s->derivative+dt/(tau+dt)*(rate-s->derivative) : 0;
    next.integral+=ki*(r-y)*dt; /* 积分绝对不能把误差也改成b*r-y。 */
    double raw=kp*(weight*r-y)+next.integral-kd*next.derivative;
    if(!isfinite(raw)) return 0;
    /* 此例不加入限幅以单独比较两自由度结构；仿真会检查输入未超过给定范围。 */
    next.previous=y; next.initialized=1; *s=next; *output=raw;
    return 1;
}
/* 对象 y''+2y'+y=u+d。临界阻尼系统在每步零阶保持输入下精确传播。 */
static void plant(double *y,double *v,double force,double dt)
{
    double z=*y-force,w=*v,decay=exp(-dt);
    *y=force+decay*(z+(w+z)*dt);
    *v=decay*(w-(w+z)*dt);
}
static double tracking(double weight,double *final)
{
    Pid s={0}; double y=0,v=0,peak=0;
    for(int k=0;k<20000;k++) {
        double u;
        assert(control(&s,1,y,.002,weight,&u));
        assert(fabs(u)<10);
        plant(&y,&v,u,.002); peak=fmax(peak,y);
    }
    *final=y;
    return peak;
}
int main(void)
{
    Pid a={0},b={0}; double ua,ub;
    assert(control(&a,1,0,.01,1,&ua));
    assert(control(&b,1,0,.01,.3,&ub));
    assert(fabs(ua-4.02)<1e-12 && fabs(ub-1.22)<1e-12);
    double end1,end2,p1=tracking(1,&end1),p2=tracking(.3,&end2);
    assert(fabs(end1-1)<.001 && fabs(end2-1)<.001 && p2<p1);
    /* 从同一物理平衡点出发：不同b需要不同积分偏置，扰动通道应相同。 */
    a=(Pid){1,1,0,1}; b=(Pid){1+4*(1-.3),1,0,1};
    double y1=1,y2=1,v1=0,v2=0,max_difference=0;
    for(int k=0;k<10000;k++) {
        double disturbance=k>=500 ? -.3 : 0;
        assert(control(&a,1,y1,.002,1,&ua));
        assert(control(&b,1,y2,.002,.3,&ub));
        plant(&y1,&v1,ua+disturbance,.002);
        plant(&y2,&v2,ub+disturbance,.002);
        max_difference=fmax(max_difference,fabs(y1-y2));
    }
    assert(max_difference<1e-10);
    assert(!control(&a,1,0,.01,2,&ua));
    printf("setpoint weighting: peak b=1 %.6f, b=0.3 %.6f; disturbance difference=%.3g\n",
           p1,p2,max_difference);
    return 0;
}

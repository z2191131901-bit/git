/* 依据已知工况变量gain查表调度PID；不是在线辨识或自适应控制。
 * 调度时重设积分贡献，补偿P、D系数变化造成的瞬时输出差。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
typedef struct { double kp,ki,kd; } Gains;
typedef struct {
    double integral,previous,derivative;
    Gains old;
    int initialized;
} Controller;
static int schedule(double gain,Gains *g)
{
    if(!g || !isfinite(gain)) return 0;
    /* 两个已整定工作点：对象增益0.25和1。区间内线性插值，区间外用端点。 */
    double t=(fmax(.25,fmin(1,gain))-.25)/.75;
    g->kp=(1-t)*8+t*2; g->ki=(1-t)*4+t*1; g->kd=(1-t)*.4+t*.1;
    return 1;
}
static int control(Controller *s,Gains g,double r,double y,double dt,double *out)
{
    if(!s || !out || !isfinite(r) || !isfinite(y) || !isfinite(dt) || dt<=0 ||
       !isfinite(g.kp) || !isfinite(g.ki) || !isfinite(g.kd) ||
       g.kp<0 || g.ki<0 || g.kd<0) return 0;
    Controller next=*s;
    double e=r-y,rate=s->initialized ? (y-s->previous)/dt : 0;
    next.derivative=s->initialized ? s->derivative+dt/(.05+dt)*(rate-s->derivative) : 0;
    /* I以输出单位保存，无需因Ki变化直接乘比例；只匹配P、D的改变。 */
    double matched=s->integral;
    if(s->initialized)
        matched+=(s->old.kp-g.kp)*e+(g.kd-s->old.kd)*next.derivative;
    double increment=g.ki*e*dt;
    next.integral=matched+increment;
    double raw=g.kp*e+next.integral-g.kd*next.derivative;
    if((raw>5 && increment>0) || (raw<-5 && increment<0)) {
        next.integral=matched; raw=g.kp*e+matched-g.kd*next.derivative;
    }
    if(!isfinite(raw) || !isfinite(next.integral)) return 0;
    next.old=g; next.previous=y; next.initialized=1;
    *s=next; *out=fmax(-5,fmin(5,raw));
    return 1;
}
static double run(int scheduled,double *final)
{
    Controller c={0}; double y=0,area=0;
    const double dt=.005,a=exp(-dt);
    for(int k=0;k<6000;k++) {
        double plant_gain=k<1000 ? 1 : .25;
        Gains g={2,1,.1}; double u;
        if(scheduled) assert(schedule(plant_gain,&g));
        assert(control(&c,g,1,y,dt,&u));
        y=a*y+(1-a)*plant_gain*u;
        if(k>=1000 && k<3000) area+=fabs(1-y)*dt;
    }
    *final=y; return area;
}
int main(void)
{
    Gains g;
    assert(schedule(.25,&g) && g.kp==8);
    assert(schedule(.625,&g) && fabs(g.kp-5)<1e-12);
    assert(schedule(2,&g) && g.kp==2);
    assert(!schedule(NAN,&g));
    /* 同一测量和误差，仅改Kp、Kd，Ki=0排除正常积分增量，验证输出匹配。 */
    Controller c={0}; double before,after;
    assert(control(&c,(Gains){2,0,.1},1,.3,.01,&before));
    assert(control(&c,(Gains){8,0,.4},1,.3,.01,&after));
    assert(fabs(after-before)<1e-12);
    double end1,end2,area1=run(0,&end1),area2=run(1,&end2);
    assert(area2<area1 && fabs(end2-1)<.001);
    printf("gain scheduling: post-change IAE fixed=%.6f scheduled=%.6f final=%.6f\n",
           area1,area2,end2);
    return 0;
}

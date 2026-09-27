/* 对称七段限加加速度轨迹：仅静止到静止，支持短距离和反向。
 * 七段jerk符号为 +,0,-,0,-,0,+；某些段持续时间可以为0。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
typedef struct { double start,goal,sign,jerk,duration[7],total; } Plan;
typedef struct { double p,v,a; } Motion;
static void advance(Motion *s,double jerk,double t)
{
    s->p+=s->v*t+s->a*t*t/2+jerk*t*t*t/6;
    s->v+=s->a*t+jerk*t*t/2;
    s->a+=jerk*t;
}
static int plan(double start,double goal,double vmax,double amax,double jmax,Plan *out)
{
    if(!out || !isfinite(start) || !isfinite(goal) || !isfinite(vmax) ||
       !isfinite(amax) || !isfinite(jmax) || vmax<=0 || amax<=0 || jmax<=0) return 0;
    double distance=fabs(goal-start);
    if(!isfinite(distance)) return 0;
    Plan p={start,goal,goal>=start?1:-1,jmax,{0},0};
    if(distance==0) { *out=p; return 1; }
    /* 先尝试达到速度上限。速度过低时，尚未达到amax就要结束加速。 */
    double tj=fmin(amax/jmax,sqrt(vmax/jmax));
    double ta=fmax(0,vmax/(jmax*tj)-tj),tv=0;
    double needed=vmax*(2*tj+ta);
    if(distance>=needed) tv=(distance-needed)/vmax;
    else {
        /* 无匀速段：判断是否还有恒加速度平台。 */
        tj=amax/jmax;
        if(distance>=2*jmax*tj*tj*tj)
            ta=fmax(0,(-3*tj+sqrt(tj*tj+4*distance/amax))/2);
        else { tj=cbrt(distance/(2*jmax)); ta=0; }
    }
    const double times[7]={tj,ta,tj,tv,tj,ta,tj};
    for(int i=0;i<7;i++) {
        if(!isfinite(times[i]) || times[i]<0) return 0;
        p.duration[i]=times[i]; p.total+=times[i];
    }
    if(!isfinite(p.total) || p.total<=0) return 0;
    *out=p; return 1;
}
static int sample(const Plan *p,double time,Motion *out)
{
    if(!p || !out || !isfinite(time) || time<0) return 0;
    if(time>=p->total) { *out=(Motion){p->goal,0,0}; return 1; }
    const int signs[7]={1,0,-1,0,-1,0,1};
    Motion s={0,0,0}; double remaining=time;
    for(int i=0;i<7 && remaining>0;i++) {
        double t=fmin(remaining,p->duration[i]);
        advance(&s,signs[i]*p->jerk,t); remaining-=t;
    }
    s.p=p->start+p->sign*s.p; s.v*=p->sign; s.a*=p->sign;
    *out=s; return 1;
}
static void check(double distance,double vmax,double amax,double jerk)
{
    Plan p; assert(plan(0,distance,vmax,amax,jerk,&p));
    Motion raw={0,0,0}; const int signs[7]={1,0,-1,0,-1,0,1};
    /* 直接积分全部段，避免末点强制置goal掩盖规划误差。 */
    for(int i=0;i<7;i++) advance(&raw,signs[i]*jerk,p.duration[i]);
    assert(fabs(raw.p-fabs(distance))<1e-9 && fabs(raw.v)<1e-9 && fabs(raw.a)<1e-9);
    double previous=0;
    for(int k=0;k<=2000;k++) {
        Motion s; assert(sample(&p,p.total*k/2000,&s));
        assert(fabs(s.v)<=vmax+1e-9 && fabs(s.a)<=amax+1e-9);
        assert(p.sign*(s.p-previous)>=-1e-9); previous=s.p;
    }
}
int main(void)
{
    check(10,2,1,2); check(2,10,1,2); check(.01,2,1,2);
    check(-2,2,1,2); check(.3,.1,2,1);
    check(.5,10,1,2); /* 恰好没有恒加速度段。 */
    Plan p; Motion s;
    assert(plan(1,1,1,1,1,&p) && p.total==0);
    assert(sample(&p,0,&s) && s.p==1 && s.v==0);
    assert(!plan(0,1,0,1,1,&p));
    assert(plan(0,10,2,1,2,&p));
    printf("S curve: PASS, 10 m total=%.6f s, jerk ramp=%.6f s, cruise=%.6f s\n",
           p.total,p.duration[0],p.duration[3]);
    return 0;
}

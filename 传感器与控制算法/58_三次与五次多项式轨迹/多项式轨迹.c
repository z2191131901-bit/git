/* 归一化时间多项式：p(s)=sum(c[i]*s^i)，s=t/T。
 * 三次指定起终点位置/速度；五次额外指定起终点加速度。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
typedef struct { double c[6],duration; } Trajectory;
typedef struct { double p,v,a,j; } Sample;
static int plan(double p0,double v0,double a0,double p1,double v1,double a1,
                double duration,int degree,Trajectory *out)
{
    if(!out || !isfinite(p0) || !isfinite(v0) || !isfinite(a0) ||
       !isfinite(p1) || !isfinite(v1) || !isfinite(a1) ||
       !isfinite(duration) || duration<=0 || (degree!=3 && degree!=5)) return 0;
    Trajectory q={{p0,v0*duration,0,0,0,0},duration};
    if(degree==3) {
        q.c[2]=3*(p1-p0)-(2*v0+v1)*duration;
        q.c[3]=-2*(p1-p0)+(v0+v1)*duration;
    } else {
        q.c[2]=a0*duration*duration/2;
        double d=p1-q.c[0]-q.c[1]-q.c[2];
        double v=v1*duration-q.c[1]-2*q.c[2];
        double a=a1*duration*duration-2*q.c[2];
        q.c[3]=10*d-4*v+a/2;
        q.c[4]=-15*d+7*v-a;
        q.c[5]=6*d-3*v+a/2;
    }
    for(int i=0;i<6;i++) if(!isfinite(q.c[i])) return 0;
    *out=q; return 1;
}
static int sample(const Trajectory *q,double time,Sample *out)
{
    if(!q || !out || !isfinite(time) || !isfinite(q->duration) ||
       q->duration<=0 || time<0 || time>q->duration) return 0;
    const double s=time/q->duration,t=q->duration;
    const double *c=q->c;
    /* 霍纳法求值；对归一化时间求导后，需要除以T、T²、T³。 */
    Sample x;
    x.p=c[0]+s*(c[1]+s*(c[2]+s*(c[3]+s*(c[4]+s*c[5]))));
    x.v=(c[1]+s*(2*c[2]+s*(3*c[3]+s*(4*c[4]+s*5*c[5]))))/t;
    x.a=(2*c[2]+s*(6*c[3]+s*(12*c[4]+s*20*c[5])))/(t*t);
    x.j=(6*c[3]+s*(24*c[4]+s*60*c[5]))/(t*t*t);
    if(!isfinite(x.p) || !isfinite(x.v) || !isfinite(x.a) || !isfinite(x.j)) return 0;
    *out=x; return 1;
}
int main(void)
{
    Trajectory cubic,quintic; Sample c,q;
    assert(plan(0,0,0,1,0,0,1,3,&cubic));
    assert(plan(0,0,0,1,0,0,1,5,&quintic));
    assert(sample(&cubic,.5,&c) && sample(&quintic,.5,&q));
    assert(fabs(c.p-.5)<1e-12 && fabs(c.v-1.5)<1e-12);
    assert(fabs(q.p-.5)<1e-12 && fabs(q.v-1.875)<1e-12);
    assert(sample(&cubic,0,&c) && sample(&quintic,0,&q));
    assert(c.a==6 && q.a==0);
    /* 非零起终点速度/加速度和反向运动，不能只验证静止到静止特例。 */
    assert(plan(2,.3,-.2,-1,-.1,.4,2,5,&quintic));
    assert(sample(&quintic,0,&q));
    assert(q.p==2 && fabs(q.v-.3)<1e-12 && fabs(q.a+.2)<1e-12);
    assert(sample(&quintic,2,&q));
    assert(fabs(q.p+1)<1e-12 && fabs(q.v+.1)<1e-12 && fabs(q.a-.4)<1e-12);
    assert(plan(2,.3,0,-1,-.1,0,2,3,&cubic));
    assert(sample(&cubic,2,&c) && fabs(c.p+1)<1e-12 && fabs(c.v+.1)<1e-12);
    assert(!sample(&cubic,3,&c));
    assert(!plan(0,0,0,1,0,0,0,5,&quintic));
    puts("polynomial: PASS, midpoint speeds cubic=1.500000 quintic=1.875000; start acceleration=6 / 0");
    return 0;
}

/*
 * 一阶离散死拍控制：y[k+1]=a*y[k]+b*u[k]。
 * 若模型精确且输出不受限，u=(r[k+1]-a*y[k])/b令下一拍输出等于目标。
 * 本例明确展示执行器限幅、模型失配和测量噪声会破坏理想条件。
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "../公共代码/PID核心.h"
typedef struct { double a,b,minimum,maximum; } Model;
typedef struct { double requested,applied; int limited; } Command;
static int deadbeat(const Model *m,double measurement,double next_reference,Command *out)
{
    if(!m || !out || !isfinite(m->a) || !isfinite(m->b) || fabs(m->b)<1e-12 ||
       !isfinite(m->minimum) || !isfinite(m->maximum) || m->minimum>=m->maximum ||
       !isfinite(measurement) || !isfinite(next_reference)) return 0;
    Command next;
    next.requested=(next_reference-m->a*measurement)/m->b;
    if(!isfinite(next.requested)) return 0;
    next.applied=control_clamp(next.requested,m->minimum,m->maximum);
    next.limited=(next.applied!=next.requested);
    *out=next; return 1;
}
typedef struct { double iae,final; int first_near,limited_steps; } Result;
static Result simulate(int use_deadbeat,double limit,double real_a,double real_b)
{
    const Model nominal={.9,.1,-limit,limit};
    const PidConfig cfg={1.5,2,0,0,-limit,limit};
    PidState pi={0};
    const double dt=.1,reference=.8;
    double y=0,area=0;
    int first=-1,saturated=0;
    for(int k=0;k<200;k++) {
        double u;
        if(use_deadbeat) {
            Command command;
            assert(deadbeat(&nominal,y,reference,&command));
            u=command.applied; saturated+=command.limited;
        } else assert(pid_update(&pi,&cfg,reference,y,0,dt,&u));
        assert(fabs(u)<=limit);
        y=real_a*y+real_b*u;
        assert(isfinite(y));
        area+=fabs(reference-y)*dt;
        /* 记录首次进入1%目标误差带，不将它称作永不离开的整定时间。 */
        if(first<0 && fabs(reference-y)<=.01*reference) first=k+1;
    }
    Result result={area,y,first,saturated}; return result;
}
int main(void)
{
    Model m={.9,.1,-100,100}; Command c;
    assert(deadbeat(&m,.2,.8,&c));
    assert(fabs(c.requested-6.2)<1e-12 && !c.limited);
    assert(fabs(.9*.2+.1*c.applied-.8)<1e-12);
    /* 负输入增益同样可解，符号由模型决定。 */
    Model negative={.9,-.1,-100,100};
    assert(deadbeat(&negative,.2,.8,&c) && fabs(c.applied+6.2)<1e-12);
    m.maximum=1; m.minimum=-1;
    assert(deadbeat(&m,0,.8,&c) && c.limited && c.applied==1);
    Command saved=c; m.b=0;
    assert(!deadbeat(&m,0,.8,&c) && c.applied==saved.applied && c.requested==saved.requested);

    Result ideal=simulate(1,100,.9,.1);
    Result limited=simulate(1,1,.9,.1);
    Result pi=simulate(0,1,.9,.1);
    Result mismatch=simulate(1,100,.85,.1); /* 放宽限幅，单独观察模型失配。 */
    assert(ideal.first_near==1 && ideal.iae<1e-12);
    assert(limited.first_near>1 && limited.limited_steps>0);
    assert(limited.iae<pi.iae && fabs(limited.final-.8)<1e-12);
    assert(fabs(pi.final-.8)<1e-6);
    /* 失配闭环：y_next=(.85-.9)*y+.8，稳态为.8/1.05。 */
    assert(fabs(mismatch.final-.8/1.05)<1e-12);
    m=(Model){.9,.1,-100,100};
    Command clean,noisy;
    assert(deadbeat(&m,.2,.8,&clean) && deadbeat(&m,.21,.8,&noisy));
    assert(fabs((noisy.applied-clean.applied)+.09)<1e-12);
    printf("deadbeat: first 1%% band step ideal=%d limited=%d PI=%d; limited steps=%d\n",
           ideal.first_near,limited.first_near,pi.first_near,limited.limited_steps);
    printf("deadbeat: 20s IAE limited=%.6f PI=%.6f; mismatched final=%.6f\n",
           limited.iae,pi.iae,mismatch.final);
    return 0;
}

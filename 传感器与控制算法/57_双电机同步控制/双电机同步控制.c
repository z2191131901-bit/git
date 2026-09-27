/* 两个速度PI外加位置差交叉反馈：快的一侧减速，慢的一侧加速。
 * 演示直线同步，不直接适用于差速转弯；转弯应跟踪期望左右位置差。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "../公共代码/PID核心.h"
typedef struct { double position,velocity; } Motor;
/* v'=-a*v+b*u+disturbance，输入保持时同时精确更新位置和速度。 */
static void plant(Motor *m,double a,double b,double u,double disturbance,double dt)
{
    double terminal=(b*u+disturbance)/a,decay=exp(-a*dt);
    m->position+=terminal*dt+(m->velocity-terminal)*(1-decay)/a;
    m->velocity=decay*m->velocity+(1-decay)*terminal;
}
static int references(double p1,double p2,double common,double gain,double limit,
                      double *r1,double *r2)
{
    if(!r1 || !r2 || !isfinite(p1) || !isfinite(p2) || !isfinite(common) ||
       !isfinite(gain) || gain<0 || !isfinite(limit) || limit<0) return 0;
    double correction=gain*(p1-p2);
    if(!isfinite(correction)) return 0;
    correction=control_clamp(correction,-limit,limit);
    double first=common-correction,second=common+correction;
    if(!isfinite(first) || !isfinite(second)) return 0;
    *r1=first; *r2=second; return 1;
}
static double run(double sync_gain,double *final_error,double *average_velocity)
{
    Motor a={0},b={0};
    PidState ca={0},cb={0};
    const PidConfig cfg={2,3,0,0,-3,3}; /* PI是PID中Kd=0的特例。 */
    double max_error=0;
    for(int k=0;k<4000;k++) {
        double r1,r2,u1,u2;
        assert(references(a.position,b.position,1,sync_gain,.3,&r1,&r2));
        assert(fabs((r1+r2)/2-1)<1e-12);
        assert(pid_update(&ca,&cfg,r1,a.velocity,0,.005,&u1));
        assert(pid_update(&cb,&cfg,r2,b.velocity,0,.005,&u2));
        /* 两台电机参数不同，第二台从第3秒开始受到额外阻力。 */
        plant(&a,2,2,u1,0,.005);
        plant(&b,2.5,1.8,u2,k>=600 ? -.4 : 0,.005);
        max_error=fmax(max_error,fabs(a.position-b.position));
    }
    *final_error=fabs(a.position-b.position);
    *average_velocity=(a.velocity+b.velocity)/2;
    return max_error;
}
int main(void)
{
    double r1,r2;
    assert(references(1,.8,1,2,.3,&r1,&r2));
    assert(fabs(r1-.7)<1e-12 && fabs(r2-1.3)<1e-12);
    assert(!references(1,.8,1,-1,.3,&r1,&r2));
    double old_error,new_error,old_speed,new_speed;
    double old_peak=run(0,&old_error,&old_speed),new_peak=run(2,&new_error,&new_speed);
    assert(old_error>.05 && new_error<.001 && new_error<old_error*.05);
    assert(new_peak<old_peak && fabs(new_speed-1)<.001);
    printf("motor sync: final position gap independent=%.6f coupled=%.6f m; peak=%.6f / %.6f\n",
           old_error,new_error,old_peak,new_peak);
    return 0;
}

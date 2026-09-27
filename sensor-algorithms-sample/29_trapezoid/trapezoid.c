/*
 * 静止到静止的梯形/三角形速度规划。输入位置m、速度m/s、加速度m/s²。
 * 规划器给出参考位置/速度/加速度，不直接驱动电机。
 */
#include "../common/pid_core.h"
#include <assert.h>
#include <stdio.h>

typedef struct
{
    double start,goal,direction,acceleration,peak_velocity;
    double acceleration_time,cruise_time,total_time;
} Trapezoid;
typedef struct { double position,velocity,acceleration; } MotionReference;

static int trajectory_plan(double start, double goal, double vmax, double amax, Trapezoid *out)
{
    if (out==NULL || !isfinite(start) || !isfinite(goal) ||
        !isfinite(vmax) || vmax<=0 || !isfinite(amax) || amax<=0) return 0;
    const double distance=fabs(goal-start);
    if (!isfinite(distance)) return 0;
    Trapezoid next={start,goal,goal>=start?1:-1,amax,0,0,0,0};
    /* 加速和减速两段总路程 v_peak²/a。距离短时峰值必须降低。 */
    next.peak_velocity=fmin(vmax,sqrt(distance)*sqrt(amax));
    next.acceleration_time=next.peak_velocity/amax;
    const double ramps_distance=next.peak_velocity*next.acceleration_time;
    next.cruise_time=next.peak_velocity>0
        ? fmax(0,(distance-ramps_distance)/next.peak_velocity) : 0;
    next.total_time=2*next.acceleration_time+next.cruise_time;
    if (!isfinite(next.total_time) || !isfinite(next.peak_velocity)) return 0;
    *out=next;
    return 1;
}
/* plan 必须来自成功的 trajectory_plan，t为从轨迹开始累计的秒数。 */
static int trajectory_sample(const Trapezoid *plan, double t, MotionReference *out)
{
    if (plan==NULL || out==NULL || !isfinite(t) || t<0) return 0;
    double s=0,v=0,a=0;
    if (t>=plan->total_time)
    {
        *out=(MotionReference){plan->goal,0,0};
        return 1;
    }
    if (t<plan->acceleration_time)
    {
        a=plan->acceleration;
        v=a*t; s=0.5*a*t*t;
    }
    else if (t<plan->acceleration_time+plan->cruise_time)
    {
        v=plan->peak_velocity;
        s=0.5*v*plan->acceleration_time+v*(t-plan->acceleration_time);
    }
    else
    {
        /* 从终点倒算剩余距离，减速末端自然到达目标且速度为零。 */
        const double remaining=plan->total_time-t;
        a=-plan->acceleration;
        v=plan->acceleration*remaining;
        s=fabs(plan->goal-plan->start)-0.5*plan->acceleration*remaining*remaining;
    }
    MotionReference next={plan->start+plan->direction*s,plan->direction*v,plan->direction*a};
    if (!isfinite(next.position) || !isfinite(next.velocity) || !isfinite(next.acceleration))
        return 0;
    *out=next;
    return 1;
}
int main(void)
{
    Trapezoid plan;
    MotionReference ref;
    int ok=trajectory_plan(0,2,1,1,&plan);
    assert(ok && fabs(plan.total_time-3)<1e-12 && plan.cruise_time==1);
    ok=trajectory_sample(&plan,0.5,&ref);
    assert(ok && fabs(ref.position-0.125)<1e-12 && ref.velocity==0.5);
    ok=trajectory_sample(&plan,2.5,&ref);
    assert(ok && fabs(ref.position-1.875)<1e-12 && ref.velocity==0.5);
    ok=trajectory_plan(2,1.75,1,1,&plan);
    assert(ok && plan.cruise_time==0 && plan.total_time==1);
    ok=trajectory_sample(&plan,0.5,&ref);
    assert(ok && fabs(ref.position-1.875)<1e-12 && ref.velocity==-0.5);
    ok=trajectory_plan(1,1,1,1,&plan);
    assert(ok);
    ok=trajectory_sample(&plan,0,&ref);
    assert(ok && ref.position==1 && ref.velocity==0);
    double saved=ref.position;
    ok=trajectory_sample(&plan,-1,&ref);
    assert(!ok && ref.position==saved);

    /* 把轨迹接到“位置P+速度PI”，用计划速度和加速度做前馈。 */
    ok=trajectory_plan(0,2,1,1,&plan);
    assert(ok);
    PidConfig cfg={4,3,0,0,-3,3};
    PidState speed_pid={0};
    double position=0,velocity=0,command=0;
    const double dt=0.01,decay=exp(-0.5*dt);
    for (int sample=0;sample<1000;++sample)
    {
        ok=trajectory_sample(&plan,sample*dt,&ref);
        assert(ok && fabs(ref.velocity)<=1+1e-12 && fabs(ref.acceleration)<=1);
        const double speed_reference=control_clamp(ref.velocity+2*(ref.position-position),-1.5,1.5);
        ok=pid_update(&speed_pid,&cfg,speed_reference,velocity,
                      ref.acceleration+0.5*ref.velocity,dt,&command);
        assert(ok && fabs(command)<=3);
        const double integral_decay=(1-decay)/0.5;
        position+=velocity*integral_decay+command/0.5*(dt-integral_decay);
        velocity=decay*velocity+command/0.5*(1-decay);
    }
    assert(fabs(position-2)<1e-4 && fabs(velocity)<1e-4);
    printf("trajectory+PI: position=%.7f velocity=%.7f\n",position,velocity);
    puts("trapezoid: PASS");
    return 0;
}

/*
 * 一阶模型逆+低通的扰动观测器（DOB教学版）。
 * 模型 y_dot=-a*y+b*u+d；d为输出变化率端的加性扰动。
 * 必须传入上一周期实际施加的u，而不是尚未限幅的期望u。
 */
#include "../公共代码/PID核心.h"
#include <assert.h>
#include <stdio.h>

typedef struct { double a_discrete,b_discrete,d_discrete,alpha; } DobConfig;
typedef struct { double previous_measurement,disturbance; int initialized; } DobState;

static int dob_design(double a, double b, double dt, double tau, DobConfig *out)
{
    if (out==NULL || !isfinite(a) || a<=0 || !isfinite(b) || b<=0 ||
        !isfinite(dt) || dt<=0 || !isfinite(tau) || tau<=0) return 0;
    /* -expm1(-a*dt) 比直接算1-exp(-a*dt)在很小dt下更准确。 */
    const double one_minus=-expm1(-a*dt);
    DobConfig next={exp(-a*dt),b*one_minus/a,one_minus/a,dt/(tau+dt)};
    if (!isfinite(next.b_discrete) || !isfinite(next.d_discrete) ||
        next.d_discrete<=0 || !isfinite(next.alpha) || next.alpha<=0) return 0;
    *out=next;
    return 1;
}
static int dob_update(DobState *state, const DobConfig *cfg, double measurement,
                       double previous_applied_input)
{
    if (state==NULL || cfg==NULL || !isfinite(measurement) ||
        !isfinite(previous_applied_input) || !isfinite(state->previous_measurement) ||
        !isfinite(state->disturbance)) return 0;
    DobState next=*state;
    if (!state->initialized)
    {
        next.previous_measurement=measurement; next.disturbance=0; next.initialized=1;
    }
    else
    {
        /* 实测输出减去“无扰动模型预测”，再除以扰动对一个周期的作用系数。 */
        const double raw=(measurement-cfg->a_discrete*state->previous_measurement
                          -cfg->b_discrete*previous_applied_input)/cfg->d_discrete;
        next.disturbance+=cfg->alpha*(raw-state->disturbance);
        next.previous_measurement=measurement;
        if (!isfinite(raw) || !isfinite(next.disturbance)) return 0;
    }
    *state=next;
    return 1;
}
static double simulate(int compensation, double *final_disturbance)
{
    DobConfig cfg;
    int ok=dob_design(2,1,0.01,0.05,&cfg);
    assert(ok);
    DobState observer={0};
    double y=0,previous_input=0;
    for (int i=0;i<2000;++i)
    {
        ok=dob_update(&observer,&cfg,y,previous_input);
        assert(ok);
        const double estimate=compensation?observer.disturbance:0;
        /* 模型前馈+P反馈-扰动补偿，参数a=2、b=1、反馈增益3。 */
        const double input=control_clamp(2*0.5+3*(0.5-y)-estimate,-3,3);
        const double disturbance=i>=500?0.4:0;
        y=cfg.a_discrete*y+cfg.b_discrete*input+cfg.d_discrete*disturbance;
        previous_input=input;
    }
    *final_disturbance=observer.disturbance;
    return y;
}
int main(void)
{
    double estimated;
    const double without=simulate(0,&estimated),with=simulate(1,&estimated);
    assert(fabs(without-0.58)<1e-8 && fabs(with-0.5)<1e-8 && fabs(estimated-0.4)<1e-8);
    DobConfig cfg;
    int ok=dob_design(2,1,0.01,0.05,&cfg);
    assert(ok);
    DobState state={0};
    ok=dob_update(&state,&cfg,0,0);
    assert(ok);
    /* 纯扰动d=.4产生的下一读数，第一次滤波估计为alpha*.4。 */
    ok=dob_update(&state,&cfg,cfg.d_discrete*0.4,0);
    assert(ok && fabs(state.disturbance-cfg.alpha*0.4)<1e-12);
    double saved=state.disturbance;
    ok=dob_update(&state,&cfg,NAN,0);
    assert(!ok && state.disturbance==saved);
    printf("DOB: baseline=%.6f compensated=%.6f disturbance=%.6f\n",without,with,estimated);
    puts("disturbance_observer: PASS");
    return 0;
}

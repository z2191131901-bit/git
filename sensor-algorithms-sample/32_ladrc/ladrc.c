/*
 * 一阶线性ADRC：LESO估计输出z1与总扰动z2，再作线性误差反馈和补偿。
 * y_dot=b0*u+f。f包含未知动力学、输入增益误差和外扰。
 * 是一阶LADRC，不是含非线性TD/fal的完整非线性ADRC。
 */
#include "../common/pid_core.h"
#include <assert.h>
#include <stdio.h>

typedef struct { double dt,b0,wc,l1,l2,minimum,maximum; } LadrcConfig;
typedef struct { double z1,z2; int initialized; } LadrcState;

static int ladrc_design(double dt, double b0, double wc, double wo,
                        double minimum, double maximum, LadrcConfig *out)
{
    if (out==NULL || !isfinite(dt) || dt<=0 || !isfinite(b0) || b0<=0 ||
        !isfinite(wc) || wc<=0 || !isfinite(wo) || wo<=0 ||
        !isfinite(minimum) || !isfinite(maximum) || minimum>=maximum) return 0;
    const double pole=exp(-wo*dt);
    if (pole<=0 || pole>=1) return 0;
    /* 对预测-校正型扩张观测器配置两个重复离散极点pole。 */
    LadrcConfig next={dt,b0,wc,1-pole*pole,(1-pole)*(1-pole)/dt,minimum,maximum};
    if (!isfinite(next.l2)) return 0;
    *out=next;
    return 1;
}
static int ladrc_update(LadrcState *state, const LadrcConfig *cfg,
                        double reference, double measurement,
                        double previous_applied_input, double *output)
{
    if (state==NULL || cfg==NULL || output==NULL || !isfinite(reference) ||
        !isfinite(measurement) || !isfinite(previous_applied_input) ||
        !isfinite(state->z1) || !isfinite(state->z2)) return 0;
    LadrcState next=*state;
    if (!state->initialized)
    {
        /* 首次没有前一段时间的数据，直接用实测初始化输出估计。 */
        next.z1=measurement; next.z2=0; next.initialized=1;
    }
    else
    {
        const double predicted=state->z1+cfg->dt*(state->z2+cfg->b0*previous_applied_input);
        const double innovation=measurement-predicted;
        next.z1=predicted+cfg->l1*innovation;
        next.z2=state->z2+cfg->l2*innovation;
    }
    const double desired_rate=cfg->wc*(reference-next.z1);
    const double raw=(desired_rate-next.z2)/cfg->b0;
    if (!isfinite(next.z1) || !isfinite(next.z2) || !isfinite(raw)) return 0;
    const double input=control_clamp(raw,cfg->minimum,cfg->maximum);
    *state=next; *output=input;
    return 1;
}
int main(void)
{
    const double dt=0.005;
    LadrcConfig cfg;
    int ok=ladrc_design(dt,1,4,20,-3,3,&cfg);
    assert(ok);
    double pole=exp(-20*dt);
    assert(fabs((2-cfg.l1-dt*cfg.l2)-2*pole)<1e-12);
    assert(fabs((1-cfg.l1)-pole*pole)<1e-12);

    LadrcState state={0};
    double y=0,previous_input=0,input=0;
    const double a=1.5,b=0.8,decay=exp(-a*dt);
    for (int i=0;i<4000;++i)
    {
        ok=ladrc_update(&state,&cfg,0.5,y,previous_input,&input);
        assert(ok && fabs(input)<=3);
        const double disturbance=i>=1000?0.4:0;
        /* 仿真真实对象与控制器b0不同，且控制器不知道这里的a。 */
        y=decay*y+(1-decay)/a*(b*input+disturbance);
        previous_input=input;
    }
    assert(fabs(y-0.5)<1e-5 && fabs(input-0.4375)<1e-5);
    assert(fabs(state.z2+0.4375)<1e-5);
    double saved=state.z2,saved_input=input;
    ok=ladrc_update(&state,&cfg,0.5,NAN,previous_input,&input);
    assert(!ok && state.z2==saved && input==saved_input);
    printf("LADRC: y=%.6f input=%.6f estimated_total_disturbance=%.6f\n",y,input,state.z2);
    puts("ladrc: PASS");
    return 0;
}

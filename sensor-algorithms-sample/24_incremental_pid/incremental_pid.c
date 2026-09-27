/*
 * 增量式 PID：计算本次输出变化量，再加到上次实际限幅输出上。
 * 固定 dt；D 对测量做二阶差分，避免目标跳变直接引入 D 冲击。
 * 本篇没有微分低通，用来明确增量公式；实际噪声下需另行处理。
 */
#include "../common/pid_core.h"
#include <assert.h>
#include <stdio.h>

typedef struct
{
    double output, previous_error;
    double previous_measurement, older_measurement;
    double dt;
    int initialized;
} IncrementalPid;

static int incremental_update(IncrementalPid *state, const PidConfig *cfg,
                               double reference, double measurement,
                               double dt, double *output)
{
    if (state==NULL || output==NULL || !pid_config_valid(cfg) ||
        cfg->derivative_tau!=0 || !isfinite(reference) || !isfinite(measurement) ||
        !isfinite(dt) || dt<=0 || !isfinite(state->output) ||
        !isfinite(state->previous_error) || !isfinite(state->previous_measurement) ||
        !isfinite(state->older_measurement) || !isfinite(state->dt))
        return 0;
    /* 这份差分公式要求周期固定，不能只把当前dt塞进历史差分。 */
    if (state->initialized && fabs(dt-state->dt)>1e-9*dt) return 0;
    IncrementalPid next=*state;
    const double error=reference-measurement;
    const double old_error=state->initialized ? state->previous_error : 0;
    const double old_y=state->initialized ? state->previous_measurement : measurement;
    const double older_y=state->initialized ? state->older_measurement : measurement;

    const double delta_p=cfg->kp*(error-old_error);
    const double delta_i=cfg->ki*dt*error;
    const double delta_d=-cfg->kd/dt*(measurement-2*old_y+older_y);
    const double candidate=state->output+delta_p+delta_i+delta_d;
    if (!isfinite(error) || !isfinite(candidate)) return 0;
    /*
     * 存的是限幅后的输出，不另藏一个持续增长的未限幅输出。
     * 这不是对所有执行器/串级场景通用的抗饱和保证。
     */
    next.output=control_clamp(candidate,cfg->minimum,cfg->maximum);
    next.previous_error=error;
    next.older_measurement=old_y;
    next.previous_measurement=measurement;
    next.dt=dt;
    next.initialized=1;
    *state=next;
    *output=next.output;
    return 1;
}
int main(void)
{
    /* 不饱和时，与相同离散形式、相同初值的位置式PID应一致。 */
    PidConfig cfg={1.5,3,0.01,0,-100,100};
    PidState position={0};
    IncrementalPid incremental={0};
    double a=0,b=0;
    for (int i=0;i<30;++i)
    {
        const double measurement=0.01*i;
        int ok=pid_update(&position,&cfg,0.5,measurement,0,0.01,&a);
        assert(ok);
        ok=incremental_update(&incremental,&cfg,0.5,measurement,0.01,&b);
        assert(ok && fabs(a-b)<1e-10);
    }
    cfg.minimum=-1; cfg.maximum=1;
    incremental=(IncrementalPid){0};
    double y=0, decay=exp(-0.01/0.5);
    for (int i=0;i<2000;++i)
    {
        int ok=incremental_update(&incremental,&cfg,0.5,y,0.01,&b);
        assert(ok && b>=-1 && b<=1);
        y=decay*y+(1-decay)*b;
    }
    assert(fabs(y-0.5)<1e-6);
    double saved=incremental.output;
    int ok=incremental_update(&incremental,&cfg,0.5,y,0.02,&b);
    assert(!ok && incremental.output==saved);
    printf("incremental PID final y=%.8f\n",y);
    puts("incremental_pid: PASS");
    return 0;
}

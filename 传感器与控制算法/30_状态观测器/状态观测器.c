/*
 * 双积分对象的离散 Luenberger 状态观测器：只测位置，估计位置与速度。
 * xhat[k+1]=A*xhat[k]+B*u[k]+L*(y[k]-C*xhat[k])。
 * 注意输出是下一时刻估计；y[k]是当前位置，不能混入y[k+1]。
 */
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

typedef struct { double dt,l_position,l_velocity; } ObserverConfig;
typedef struct { double position,velocity; } ObserverState;

static int observer_design(double dt, double pole1, double pole2, ObserverConfig *out)
{
    if (out==NULL || !isfinite(dt) || dt<=0 || !isfinite(pole1) ||
        !isfinite(pole2) || pole1<=0 || pole1>=1 || pole2<=0 || pole2>=1)
        return 0;
    /* 令 A-LC 的特征多项式等于 (lambda-pole1)*(lambda-pole2)。 */
    ObserverConfig next={dt,2-pole1-pole2,(1-pole1)*(1-pole2)/dt};
    if (!isfinite(next.l_velocity)) return 0;
    *out=next;
    return 1;
}
static int observer_update(ObserverState *state, const ObserverConfig *cfg,
                            double measured_position, int valid, double applied_acceleration)
{
    if (state==NULL || cfg==NULL || !isfinite(state->position) ||
        !isfinite(state->velocity) || !isfinite(applied_acceleration) ||
        (valid && !isfinite(measured_position))) return 0;
    /* cfg 必须来自成功的 observer_design。缺测时只做模型预测。 */
    const double residual=valid ? measured_position-state->position : 0;
    ObserverState next={
        state->position+cfg->dt*state->velocity+0.5*cfg->dt*cfg->dt*applied_acceleration
            +cfg->l_position*residual,
        state->velocity+cfg->dt*applied_acceleration+cfg->l_velocity*residual
    };
    if (!isfinite(next.position) || !isfinite(next.velocity)) return 0;
    *state=next;
    return 1;
}
int main(void)
{
    ObserverConfig cfg;
    int ok=observer_design(0.02,0.8,0.85,&cfg);
    assert(ok && fabs(cfg.l_position-0.35)<1e-12 && fabs(cfg.l_velocity-1.5)<1e-12);
    /* trace/determinant 分别对应两个指定极点的和与积。 */
    assert(fabs((2-cfg.l_position)-(0.8+0.85))<1e-12);
    assert(fabs((1-cfg.l_position+cfg.dt*cfg.l_velocity)-0.8*0.85)<1e-12);

    ObserverState estimate={-2,1};
    double position=1,velocity=-0.2;
    for (int i=0;i<300;++i)
    {
        double input=0.3*sin(i*0.02);
        ok=observer_update(&estimate,&cfg,position,1,input);
        assert(ok);
        position+=cfg.dt*velocity+0.5*cfg.dt*cfg.dt*input;
        velocity+=cfg.dt*input;
    }
    assert(fabs(estimate.position-position)<1e-9 && fabs(estimate.velocity-velocity)<1e-9);
    estimate=(ObserverState){1,2};
    ok=observer_update(&estimate,&cfg,NAN,0,3);
    assert(ok && fabs(estimate.position-1.0406)<1e-12 && fabs(estimate.velocity-2.06)<1e-12);
    double saved=estimate.position;
    ok=observer_update(&estimate,&cfg,NAN,1,0);
    assert(!ok && estimate.position==saved);
    ok=observer_design(0.02,1.1,0.8,&cfg);
    assert(!ok);
    puts("state_observer: PASS");
    return 0;
}

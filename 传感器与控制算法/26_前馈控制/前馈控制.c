/*
 * 模型前馈+P反馈：对象 y_dot=-a*y+b*u。
 * 想跟随 reference，前馈 u_ff=(reference_dot+a*reference)/b。
 * P反馈修正跟踪误差；本例刻意不加I，展示前馈的作用而非掩盖模型差异。
 */
#include "../公共代码/PID核心.h"
#include <assert.h>
#include <stdio.h>

static double simulate(int enabled, double *accumulated_error)
{
    const double a=2,b=2,dt=0.01;
    PidConfig cfg={1,0,0,0,-1,1};
    PidState state={0};
    double y=0,u=0,total=0;
    for (int sample=0;sample<1500;++sample)
    {
        const double t=sample*dt;
        const double reference=t<5 ? 0.1*t : 0.5;
        const double reference_rate=t<5 ? 0.1 : 0;
        const double feedforward=enabled ? (reference_rate+a*reference)/b : 0;
        int ok=pid_update(&state,&cfg,reference,y,feedforward,dt,&u);
        assert(ok && fabs(u)<=1);
        total+=fabs(reference-y)*dt;
        /* 零阶保持的精确离散模型。 */
        const double decay=exp(-a*dt);
        y=decay*y+(b/a)*(1-decay)*u;
    }
    *accumulated_error=total;
    return y;
}
int main(void)
{
    double error_without,error_with;
    double without=simulate(0,&error_without);
    double with=simulate(1,&error_with);
    assert(fabs(without-0.25)<1e-6 && fabs(with-0.5)<1e-6);
    assert(error_with<0.1*error_without);
    printf("P only final=%.6f; feedforward+P final=%.6f\n",without,with);
    printf("integrated absolute error: without=%.6f with=%.6f\n",error_without,error_with);
    puts("feedforward: PASS");
    return 0;
}

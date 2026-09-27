/* 位置式 PID 示例：手算、抗积分饱和、无微分设定冲击、闭环跟踪。 */
#include "../公共代码/PID核心.h"
#include <assert.h>
#include <stdio.h>

int main(void)
{
    PidConfig cfg={2,1,0,0,-100,100};
    PidState state={0};
    double u=0;
    /* error=1，P=2，I=1*1*0.1=0.1，首次D=0，所以输出2.1。 */
    int ok=pid_update(&state,&cfg,1,0,0,0.1,&u);
    assert(ok && fabs(u-2.1)<1e-12);
    printf("hand calculation output=%.6f\n",u);

    cfg=(PidConfig){2,4,0.05,0.02,-1,1};
    state=(PidState){0};
    for (int i=0;i<1000;++i)
    {
        ok=pid_update(&state,&cfg,10,0,0,0.01,&u);
        assert(ok && u==1 && state.integral==0);
    }
    ok=pid_update(&state,&cfg,0,0,0,0.01,&u);
    assert(ok && u==0); /* 不会因长期不可达目标积累很大积分。 */

    PidConfig derivative_only={0,0,1,0,-100,100};
    PidState d_state={0};
    ok=pid_update(&d_state,&derivative_only,0,0,0,0.01,&u);
    assert(ok);
    ok=pid_update(&d_state,&derivative_only,10,0,0,0.01,&u);
    assert(ok && u==0); /* 目标跳变、测量不变，D不产生冲击。 */
    ok=pid_update(&d_state,&derivative_only,10,0.1,0,0.01,&u);
    assert(ok && fabs(u+10)<1e-12); /* 测量上升产生负向阻尼。 */

    /* 单独核对微分低通：tau=.09、dt=.01时alpha=.1。 */
    derivative_only.derivative_tau=0.09;
    d_state=(PidState){0};
    ok=pid_update(&d_state,&derivative_only,0,0,0,0.01,&u);
    assert(ok);
    ok=pid_update(&d_state,&derivative_only,0,1,0,0.01,&u);
    assert(ok && fabs(d_state.derivative-10)<1e-12 && fabs(u+10)<1e-12);

    /* 仿真对象 y_dot=(-y+u)/0.5；零阶保持输入下精确离散。 */
    state=(PidState){0};
    double y=0;
    const double dt=0.01, decay=exp(-dt/0.5);
    for (int i=0;i<2000;++i)
    {
        ok=pid_update(&state,&cfg,0.5,y,0,dt,&u);
        assert(ok && u>=-1 && u<=1);
        y=decay*y+(1-decay)*u;
    }
    assert(fabs(y-0.5)<1e-6);
    double saved_integral=state.integral, saved_u=u;
    ok=pid_update(&state,&cfg,0.5,y,0,0,&u);
    assert(!ok && state.integral==saved_integral && u==saved_u);
    printf("PID final y=%.8f reference=0.5\n",y);
    puts("pid: PASS");
    return 0;
}

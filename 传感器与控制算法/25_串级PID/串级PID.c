/*
 * 串级控制：位置外环产生速度目标，速度内环产生驱动力/加速度命令。
 * 外环为P（PID的Ki/Kd为0），内环为PI；先调内环，再调外环。
 * 外环20Hz，内环100Hz；外环速度目标在两次计算间保持。
 */
#include "../公共代码/PID核心.h"
#include <assert.h>
#include <stdio.h>

int main(void)
{
    PidConfig position_cfg={2,0,0,0,-1,1}; /* 输出单位 m/s。 */
    PidConfig velocity_cfg={4,3,0,0,-3,3}; /* 输出按单位质量归一为 m/s²。 */
    PidState position_pid={0},velocity_pid={0};
    double position=0,velocity=0,velocity_reference=0,command=0;
    double peak_position=0;
    const double dt=0.01;
    for (int sample=0;sample<2000;++sample)
    {
        if (sample%5==0)
        {
            /* 外环接收位置误差，输出是内环的目标，不能直接发给电机。 */
            int ok=pid_update(&position_pid,&position_cfg,1,position,0,5*dt,&velocity_reference);
            assert(ok && fabs(velocity_reference)<=1);
        }
        int ok=pid_update(&velocity_pid,&velocity_cfg,velocity_reference,velocity,0,dt,&command);
        assert(ok && fabs(command)<=3);

        /*
         * 仿真对象：v_dot=command-0.5*v，p_dot=v。
         * 对一个周期内不变的command使用精确离散，避免把仿真积分误差当控制效果。
         */
        const double decay=exp(-0.5*dt);
        const double integral_decay=(1-decay)/0.5;
        position+=velocity*integral_decay+command/0.5*(dt-integral_decay);
        velocity=decay*velocity+command/0.5*(1-decay);
        if (position>peak_position) peak_position=position;
    }
    assert(fabs(position-1)<1e-5 && fabs(velocity)<1e-5);
    assert(peak_position<1.5);
    printf("cascade final position=%.8f velocity=%.8f peak=%.5f\n",
           position,velocity,peak_position);
    puts("cascade_pid: PASS");
    return 0;
}

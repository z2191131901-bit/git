/*
 * 编码器测速：M 法（固定时间数脉冲）与 T 法（测脉冲间隔）。
 * counts_per_rev 指最终解码后每一轮转一圈的计数，必须包含倍频和减速比。
 * 已经由驱动处理计数器回绕；此处接收带符号的安全计数增量。
 */
#include "../common/math3d.h"
#include <assert.h>
#include <stdint.h>
#include <stdio.h>

static int speed_from_count(int64_t delta_count, double dt, double counts_per_rev,
                            double wheel_radius, double *speed)
{
    if (speed==NULL || !isfinite(dt) || dt<=0 ||
        !isfinite(counts_per_rev) || counts_per_rev<=0 ||
        !isfinite(wheel_radius) || wheel_radius<=0) return 0;
    /* 每计数对应路程=2*pi*半径/每圈计数，乘增量再除时间。 */
    const double next=((double)delta_count/counts_per_rev)*2*SENSOR_PI*wheel_radius/dt;
    if (!isfinite(next)) return 0;
    *speed=next;
    return 1;
}
static int speed_from_period(double period, int direction, double age, double timeout,
                             double counts_per_rev, double radius, double *speed)
{
    if (speed==NULL || !isfinite(period) || period<=0 ||
        (direction!=1 && direction!=-1) || !isfinite(age) || age<0 ||
        !isfinite(timeout) || timeout<=0 || !isfinite(counts_per_rev) ||
        counts_per_rev<=0 || !isfinite(radius) || radius<=0) return 0;
    /* 超过约定时间仍无新脉冲，应用层将速度置零；这是判断策略，不是精确测量。 */
    const double next=age>timeout ? 0 : direction*2*SENSOR_PI*radius/(counts_per_rev*period);
    if (!isfinite(next)) return 0;
    *speed=next;
    return 1;
}
int main(void)
{
    double speed=0;
    int ok=speed_from_count(100,0.1,1000,0.05,&speed);
    assert(ok && fabs(speed-SENSOR_PI/10)<1e-12);
    printf("count method speed=%.6f m/s\n",speed);
    ok=speed_from_count(-100,0.1,1000,0.05,&speed);
    assert(ok && speed<0);
    ok=speed_from_period(0.001,1,0.002,0.1,1000,0.05,&speed);
    assert(ok && fabs(speed-SENSOR_PI/10)<1e-12);
    ok=speed_from_period(0.001,1,0.2,0.1,1000,0.05,&speed);
    assert(ok && speed==0);
    ok=speed_from_count(1,0,1000,0.05,&speed);
    assert(!ok && speed==0);
    puts("encoder_speed: PASS");
    return 0;
}

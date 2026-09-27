/* 四元数转 Z-Y-X 欧拉角：用于显示，不建议转回三个角独立积分。 */
#include "../公共代码/三维数学.h"
#include <assert.h>
#include <stdio.h>

typedef struct { double roll, pitch, yaw; } Euler;

static int quaternion_to_euler(Quat q, Euler *output)
{
    if (output == NULL || !quat_normalize(&q)) return 0;
    const double sin_pitch = clamp_value(2*(q.w*q.y-q.z*q.x),-1,1);
    Euler next;
    next.pitch = asin(sin_pitch); /* clamp 避免舍入后超出 asin 定义域。 */
    if (fabs(sin_pitch) > 1-1e-10)
    {
        /*
         * 俯仰约 ±90 度时，roll/yaw 无法唯一分离。
         * 这里固定 roll=0，选取能重建同一方向的 yaw。
         */
        next.roll = 0;
        const double r01 = 2*(q.x*q.y-q.w*q.z);
        const double r11 = 1-2*(q.x*q.x+q.z*q.z);
        next.yaw = atan2(-r01,r11);
    }
    else
    {
        next.roll = atan2(2*(q.w*q.x+q.y*q.z),1-2*(q.x*q.x+q.y*q.y));
        next.yaw = atan2(2*(q.w*q.z+q.x*q.y),1-2*(q.y*q.y+q.z*q.z));
    }
    *output=next;
    return 1;
}
int main(void)
{
    Euler e;
    Quat q=quat_from_euler(0.2,-0.3,0.7);
    int ok=quaternion_to_euler(q,&e);
    assert(ok && fabs(e.roll-0.2)<1e-12 && fabs(e.pitch+0.3)<1e-12 && fabs(e.yaw-0.7)<1e-12);

    for (int sign=-1;sign<=1;sign+=2)
    {
        q=quat_from_euler(0.4,sign*SENSOR_PI/2,0.8);
        ok=quaternion_to_euler(q,&e);
        assert(ok && e.roll==0);
        Quat rebuilt=quat_from_euler(e.roll,e.pitch,e.yaw);
        /* q 和 -q 表示同一姿态，所以检查点积绝对值是否接近 1。 */
        double dot=q.w*rebuilt.w+q.x*rebuilt.x+q.y*rebuilt.y+q.z*rebuilt.z;
        assert(fabs(fabs(dot)-1)<1e-10);
    }
    Quat bad={0,0,0,0};
    ok=quaternion_to_euler(bad,&e);
    assert(!ok);
    puts("euler_angles: PASS");
    return 0;
}

/*
 * 经典 Madgwick 的 IMU（六轴）版本：陀螺仪+加速度计。
 * 六轴不能确定绝对航向。本篇不是九轴 MARG，也不是 Fusion 的修订算法。
 */
#include "../common/math3d.h"
#include <assert.h>
#include <stdio.h>

static int madgwick_update(Quat *state, Vec3 gyro, Vec3 accel, int accel_valid,
                           double beta, double dt)
{
    if (state==NULL || !vec_finite(gyro) || !isfinite(beta) || beta<0 ||
        !isfinite(dt) || dt<=0) return 0;
    Quat q=*state;
    if (!quat_normalize(&q)) return 0;
    Quat angular={0,gyro.x,gyro.y,gyro.z};
    Quat derivative=quat_multiply(q,angular);
    derivative.w*=0.5; derivative.x*=0.5; derivative.y*=0.5; derivative.z*=0.5;

    if (accel_valid)
    {
        Vec3 a;
        if (!vec_unit(accel,&a)) return 0;
        /*
         * f 是“根据姿态预测的向上方向 - 加速度测得的向上方向”。
         * 与本库机体到世界的姿态、静止 +Z 比力约定相配。
         */
        const double f1=2*(q.x*q.z-q.w*q.y)-a.x;
        const double f2=2*(q.w*q.x+q.y*q.z)-a.y;
        const double f3=1-2*(q.x*q.x+q.y*q.y)-a.z;
        /* gradient = J^T*f：四个分量对应对 w、x、y、z 的偏导。 */
        Quat gradient={
            -2*q.y*f1+2*q.x*f2,
             2*q.z*f1+2*q.w*f2-4*q.x*f3,
            -2*q.w*f1+2*q.z*f2-4*q.y*f3,
             2*q.x*f1+2*q.y*f2
        };
        const double length=hypot(hypot(gradient.w,gradient.x),hypot(gradient.y,gradient.z));
        /* 误差已经为零时不能除以零，此时只保留陀螺仪导数。 */
        if (length>1e-12)
        {
            derivative.w-=beta*gradient.w/length;
            derivative.x-=beta*gradient.x/length;
            derivative.y-=beta*gradient.y/length;
            derivative.z-=beta*gradient.z/length;
        }
    }
    Quat next={q.w+derivative.w*dt,q.x+derivative.x*dt,
               q.y+derivative.y*dt,q.z+derivative.z*dt};
    if (!quat_normalize(&next)) return 0;
    *state=next;
    return 1;
}
int main(void)
{
    Quat q=quat_from_euler(0.3,0,0);
    Vec3 zero={0,0,0}, up={0,0,1};
    int ok;
    for (int i=0;i<4000;++i)
    {
        ok=madgwick_update(&q,zero,up,1,0.1,0.002);
        assert(ok);
    }
    /* 固定归一化梯度步长会在最优点附近微振荡，容差随 beta*dt 选择。 */
    assert(fabs(q.x)<0.001);
    q=(Quat){1,0,0,0};
    ok=madgwick_update(&q,zero,up,1,0.1,0.01);
    assert(ok && q.w==1 && q.x==0);

    Vec3 yaw_rate={0,0,SENSOR_PI/2};
    for (int i=0;i<1000;++i)
    {
        ok=madgwick_update(&q,yaw_rate,zero,0,0.1,0.001);
        assert(ok);
    }
    assert(fabs(q.z-sqrt(0.5))<1e-6);
    Quat saved=q;
    ok=madgwick_update(&q,zero,zero,1,0.1,0.01);
    assert(!ok && q.w==saved.w);
    puts("madgwick: PASS");
    return 0;
}

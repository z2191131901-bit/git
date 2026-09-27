/*
 * Mahony 向量反馈姿态融合：加速度提供倾角参考，可选已知世界磁场方向。
 * q 为机体到世界；accel 输入为比力，水平静止时指向 +Z。
 * 磁场参考与测量必须在一致轴约定下，磁力计需先完成校准。
 */
#include "../公共代码/三维数学.h"
#include <assert.h>
#include <stdio.h>

typedef struct
{
    Quat q;
    Vec3 integral; /* 积分校正角速度，是负零偏估计的意义，单位 rad/s。 */
} Mahony;

static int mahony_update(Mahony *state, Vec3 gyro, Vec3 accel, int accel_valid,
                         Vec3 mag, Vec3 magnetic_world, int mag_valid,
                         double kp, double ki, double dt)
{
    if (state==NULL || !vec_finite(gyro) || !isfinite(kp) || kp<0 ||
        !isfinite(ki) || ki<0 || !isfinite(dt) || dt<=0) return 0;
    Mahony next=*state; /* 在本地计算，任何失败都不破坏原状态。 */
    if (!quat_normalize(&next.q) || !vec_finite(next.integral)) return 0;
    Vec3 error={0,0,0};
    if (accel_valid)
    {
        Vec3 measured;
        if (!vec_unit(accel,&measured)) return 0;
        Vec3 up={0,0,1};
        Vec3 predicted=quat_rotate(quat_conjugate(next.q),up);
        /* 测量×预测，顺序决定反馈是收敛还是发散。 */
        error=vec_add(error,vec_cross(measured,predicted));
    }
    if (mag_valid)
    {
        Vec3 measured, reference;
        if (!vec_unit(mag,&measured) || !vec_unit(magnetic_world,&reference)) return 0;
        Vec3 up={0,0,1};
        /* 磁场若与竖直近乎平行，不能可靠提供航向参考。 */
        if (vec_norm(vec_cross(reference,up))<1e-3) return 0;
        Vec3 predicted=quat_rotate(quat_conjugate(next.q),reference);
        error=vec_add(error,vec_cross(measured,predicted));
    }
    if (ki==0) next.integral=(Vec3){0,0,0};
    else if (accel_valid || mag_valid)
    {
        Vec3 change=vec_scale(error,ki*dt);
        if (!vec_finite(change)) return 0;
        next.integral=vec_add(next.integral,change);
        /* 教学防积分饱和：每轴限制 ±0.5 rad/s，工程中按传感器配置。 */
        next.integral.x=clamp_value(next.integral.x,-0.5,0.5);
        next.integral.y=clamp_value(next.integral.y,-0.5,0.5);
        next.integral.z=clamp_value(next.integral.z,-0.5,0.5);
    }
    Vec3 corrected=vec_add(gyro,vec_add(vec_scale(error,kp),next.integral));
    if (!quat_step(&next.q,corrected,dt)) return 0;
    *state=next;
    return 1;
}
int main(void)
{
    Vec3 zero={0,0,0}, up={0,0,1}, north={1,0,0};
    Mahony state={quat_from_euler(0.3,0,0),{0,0,0}};
    int ok;
    for (int i=0;i<2000;++i)
    {
        ok=mahony_update(&state,zero,up,1,zero,zero,0,2,0,0.005);
        assert(ok);
    }
    assert(fabs(state.q.x)<1e-7); /* 零陀螺输入也必须能完成倾角校正。 */

    state=(Mahony){quat_from_euler(0,0,0.4),{0,0,0}};
    for (int i=0;i<2000;++i)
    {
        ok=mahony_update(&state,zero,up,1,north,north,1,2,0,0.005);
        assert(ok);
    }
    assert(fabs(state.q.z)<1e-7); /* 已知非共线磁场可校正航向。 */

    state=(Mahony){{1,0,0,0},{0,0,0}};
    Vec3 biased={0.05,0,0};
    for (int i=0;i<12000;++i)
    {
        ok=mahony_update(&state,biased,up,1,zero,zero,0,2,0.5,0.005);
        assert(ok);
    }
    assert(fabs(state.integral.x+0.05)<1e-6 && fabs(state.q.x)<1e-6);
    Quat saved=state.q;
    ok=mahony_update(&state,zero,zero,1,zero,zero,0,2,0.5,0.005);
    assert(!ok && saved.w==state.q.w);
    printf("estimated correction x=%.6f rad/s\n",state.integral.x);
    puts("mahony: PASS");
    return 0;
}

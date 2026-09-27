/*
 * 简化局部三维捷联惯导：陀螺仪 -> 四元数 -> 比力旋转 -> 速度/位置积分。
 * q 为机体到世界，世界 Z 向上；加速度计静止水平输出 +g。
 * 输入区间平均角速度和比力；忽略地球自转、曲率、科里奥利效应。
 * 这是短时、小范围的惯性递推教学实现，没有误差状态滤波器。
 */
#include <assert.h>
#include <stdio.h>
#include "../公共代码/三维数学.h"
typedef struct { Quat q; Vec3 v, p; } Ins;
static int ins_step(Ins *state, Vec3 gyro, Vec3 force, Vec3 gyro_bias,
                    Vec3 accel_bias, double dt)
{
    if (!state || !vec_finite(state->p) || !vec_finite(state->v) ||
        !vec_finite(gyro) || !vec_finite(force) || !vec_finite(gyro_bias) ||
        !vec_finite(accel_bias) || !isfinite(dt) || dt<=0) return 0;
    Ins next=*state;
    Vec3 rate=vec_add(gyro,vec_scale(gyro_bias,-1));
    Vec3 corrected=vec_add(force,vec_scale(accel_bias,-1));
    Quat middle=state->q;
    /* 用区间中间的姿态旋转比力，降低边转动边加速的积分误差。 */
    if (!quat_step(&middle,rate,dt/2) || !quat_step(&next.q,rate,dt)) return 0;
    Vec3 acceleration=quat_rotate(middle,corrected);
    acceleration.z-=9.80665; /* a_world = R(q)*f_body + [0,0,-g]。 */
    /* 先使用旧速度算位置，再更新速度；避免误把 v_next 当 v_old。 */
    next.p=vec_add(state->p,vec_add(vec_scale(state->v,dt),
                                   vec_scale(acceleration,0.5*dt*dt)));
    next.v=vec_add(state->v,vec_scale(acceleration,dt));
    if (!vec_finite(next.p) || !vec_finite(next.v)) return 0;
    *state=next;
    return 1;
}
/* 只有外部可靠确认静止时才调用。直接置零是演示，不是卡尔曼 ZUPT。 */
static void confirmed_stop(Ins *state) { state->v=(Vec3){0,0,0}; }

int main(void)
{
    const Vec3 zero={0,0,0}, gravity={0,0,9.80665};
    Ins level={{1,0,0,0},{0,0,0},{0,0,0}};
    for (int i=0;i<1000;i++)
        assert(ins_step(&level,zero,gravity,zero,zero,.01));
    assert(vec_norm(level.p)<1e-12 && vec_norm(level.v)<1e-12);
    /* 倾斜静止：把世界 +g 旋转到机体，结果也应该保持不动。 */
    Ins tilted={quat_from_euler(.3,-.2,.7),{0,0,0},{0,0,0}};
    Vec3 body_g=quat_rotate(quat_conjugate(tilted.q),gravity);
    for (int i=0;i<100;i++)
        assert(ins_step(&tilted,zero,body_g,zero,zero,.01));
    assert(vec_norm(tilted.p)<1e-12);
    Ins accelerating={{1,0,0,0},{0,0,0},{0,0,0}};
    for (int i=0;i<100;i++)
        assert(ins_step(&accelerating,zero,(Vec3){1,0,9.80665},zero,zero,.01));
    assert(fabs(accelerating.p.x-.5)<1e-12);
    assert(fabs(accelerating.v.x-1)<1e-12);
    Ins rotating={{1,0,0,0},{0,0,0},{0,0,0}};
    for (int i=0;i<100;i++)
        assert(ins_step(&rotating,(Vec3){0,0,1},gravity,zero,zero,.01));
    assert(fabs(rotating.q.w-cos(.5))<1e-12 && vec_norm(rotating.p)<1e-12);
    Ins drifting={{1,0,0,0},{0,0,0},{0,0,0}};
    for (int i=0;i<1000;i++)
        assert(ins_step(&drifting,zero,(Vec3){.01,0,9.80665},zero,zero,.01));
    assert(fabs(drifting.p.x-.5)<1e-10);
    const double drift=drifting.p.x;
    confirmed_stop(&drifting);
    assert(vec_norm(drifting.v)==0 && drifting.p.x==drift);
    Ins saved=drifting;
    assert(!ins_step(&drifting,zero,gravity,zero,zero,0));
    assert(drifting.p.x==saved.p.x);
    printf("strapdown: 1s acceleration p=%.6f v=%.6f; 10s accel-bias drift=%.6f m\n",
           accelerating.p.x,accelerating.v.x,drift);
    return 0;
}

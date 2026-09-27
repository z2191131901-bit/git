/*
 * 坐标变换与去重力：a_world = R(q)*specific_force_body + gravity_world。
 * 世界 Z 向上，因此 gravity_world=[0,0,-9.80665]。
 * 输入必须为仍含重力相关静止分量的加速度计比力，单位 m/s^2。
 */
#include "../common/math3d.h"
#include <assert.h>
#include <stdio.h>

static int linear_acceleration(Quat q, Vec3 specific_force, double gravity, Vec3 *out)
{
    if (out==NULL || !vec_finite(specific_force) ||
        !isfinite(gravity) || gravity<=0 || !quat_normalize(&q)) return 0;
    Vec3 world=quat_rotate(q,specific_force);
    world.z-=gravity; /* 静止时旋转后为 +g，此处抵消为 0。 */
    if (!vec_finite(world)) return 0;
    *out=world;
    return 1;
}
int main(void)
{
    const double g=9.80665;
    Quat q=quat_from_euler(0.5,-0.4,0.7);
    Vec3 up_force={0,0,g};
    Vec3 body=quat_rotate(quat_conjugate(q),up_force), acceleration;
    int ok=linear_acceleration(q,body,g,&acceleration);
    assert(ok && vec_norm(acceleration)<1e-12);

    /* 已知世界向前加速 2 m/s^2，合成传感器比力应为 [2,0,g] 旋转到机体。 */
    body=quat_rotate(quat_conjugate(q),(Vec3){2,0,g});
    ok=linear_acceleration(q,body,g,&acceleration);
    assert(ok && fabs(acceleration.x-2)<1e-12 && fabs(acceleration.z)<1e-12);

    /* 自由落体比力约为零，世界加速度应向下 g，不是零。 */
    ok=linear_acceleration(q,(Vec3){0,0,0},g,&acceleration);
    assert(ok && fabs(acceleration.z+g)<1e-12);
    ok=linear_acceleration((Quat){0,0,0,0},body,g,&acceleration);
    assert(!ok);
    puts("gravity_removal: PASS");
    return 0;
}

/* 四元数乘法、共轭、归一化的可运行示例；具体公式在公共头文件。 */
#include "../common/math3d.h"
#include <assert.h>
#include <stdio.h>

int main(void)
{
    /* [2,0,0,0] 不是单位四元数，先把长度变成 1。 */
    Quat q = {2,0,0,0};
    int ok = quat_normalize(&q);
    assert(ok && q.w == 1);

    /* 绕 Z 转 90 度：[cos45,0,0,sin45]。 */
    q = quat_from_euler(0,0,SENSOR_PI/2);
    Quat inverse = quat_conjugate(q);
    Quat identity = quat_multiply(q,inverse);
    assert(fabs(identity.w-1)<1e-12 && fabs(identity.z)<1e-12);

    /* X 方向向量旋转后成为 Y 方向。这个结果比仅检查长度更有意义。 */
    Vec3 x_axis = {1,0,0};
    Vec3 y_axis = quat_rotate(q,x_axis);
    assert(fabs(y_axis.x)<1e-12 && fabs(y_axis.y-1)<1e-12);

    /* 先绕 X 再绕 Z，对应 qz*qx；交换顺序一般得不到同一结果。 */
    Quat qx = quat_from_euler(SENSOR_PI/2,0,0);
    Vec3 first = quat_rotate(quat_multiply(q,qx),x_axis);
    Vec3 second = quat_rotate(quat_multiply(qx,q),x_axis);
    assert(fabs(first.y-1)<1e-12 && fabs(second.z-1)<1e-12);

    Quat zero = {0,0,0,0};
    ok = quat_normalize(&zero);
    assert(!ok);
    printf("X rotated by yaw 90: (%.3f, %.3f, %.3f)\n",y_axis.x,y_axis.y,y_axis.z);
    puts("quaternion_basics: PASS");
    return 0;
}

/* 增量四元数更新：采样间隔内角速度固定时的旋转指数映射。 */
#include "../公共代码/三维数学.h"
#include <assert.h>
#include <stdio.h>

int main(void)
{
    Quat q = {1,0,0,0};
    Vec3 rate = {0,0,SENSOR_PI/2}; /* 90 度/秒，注意输入实际为 rad/s。 */
    int ok = quat_step(&q,rate,1.0);
    assert(ok && fabs(q.w-sqrt(0.5))<1e-12 && fabs(q.z-sqrt(0.5))<1e-12);

    /*
     * 同一恒定角速度分成 100 小步，理论结果应与一步一致。
     * 这是恒定角速度测试；并不意味着任意快速变向都能无误差恢复。
     */
    Quat small_steps = {1,0,0,0};
    for (int i=0;i<100;++i)
    {
        ok = quat_step(&small_steps,rate,0.01);
        assert(ok);
    }
    assert(fabs(small_steps.w-q.w)<1e-12 && fabs(small_steps.z-q.z)<1e-12);
    Vec3 zero = {0,0,0};
    ok = quat_step(&q,zero,0.01);
    assert(ok && fabs(q.z-sqrt(0.5))<1e-12);
    const Quat saved = q;
    ok = quat_step(&q,rate,-0.1);
    assert(!ok && q.w==saved.w && q.z==saved.z);
    printf("q=(%.6f,%.6f,%.6f,%.6f)\n",q.w,q.x,q.y,q.z);
    puts("quaternion_increment: PASS");
    return 0;
}

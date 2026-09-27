/*
 * 加速度计倾角：仅在静止或平动加速度可忽略时有意义。
 * 右手机体坐标系 X 前、Y 左、Z 上，水平静止读数为 [0,0,+g]。
 * 姿态使用 Z-Y-X 欧拉角约定；输出为弧度，不估计航向。
 */
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

/* 输出结构体把两个相关结果放在一起。 */
typedef struct
{
    float roll;  /* 横滚：绕 X 轴。 */
    float pitch; /* 俯仰：绕 Y 轴，主值范围 [-pi/2, pi/2]。 */
} Tilt;

static int accel_tilt(float ax, float ay, float az, Tilt *output)
{
    if (output == NULL || !isfinite(ax) || !isfinite(ay) || !isfinite(az))
        return 0;

    /*
     * hypot(a,b) = sqrt(a*a+b*b)，内部处理对大/小数更稳健。
     * yz_length 是加速度在 YZ 平面上的投影长度。
     */
    const double yz_length = hypot((double)ay, (double)az);
    const double length = hypot((double)ax, yz_length);

    /*
     * 零向量没有方向；接近 X 轴时横滚不可可靠确定。
     * 用比例判断奇异区域，避免依赖输入使用 g 还是 m/s^2。
     */
    if (length == 0 || yz_length / length < 1e-6)
        return 0;

    Tilt next;
    /* atan2 同时使用两个参数的符号区分象限，不能简单换成 atan(ay/az)。 */
    next.roll = (float)atan2((double)ay, (double)az);
    next.pitch = (float)atan2(-(double)ax, yz_length);
    *output = next;
    return 1;
}

int main(void)
{
    const float pi = 3.14159265358979323846f;
    Tilt tilt = {0, 0};
    int ok = accel_tilt(0, 0, 1, &tilt);
    assert(ok && tilt.roll == 0 && tilt.pitch == 0);

    /* 横滚 +30 度、俯仰 0 度，静止加速度方向为 [0,sin30,cos30]。 */
    ok = accel_tilt(0, 0.5f, 0.8660254f, &tilt);
    assert(ok && fabsf(tilt.roll - pi / 6) < 1e-5f);
    printf("roll=%.3f deg\n", tilt.roll * 180 / pi);

    /* 俯仰 +30 度，对应 X 分量为负。 */
    ok = accel_tilt(-0.5f, 0, 0.8660254f, &tilt);
    assert(ok && fabsf(tilt.pitch - pi / 6) < 1e-5f);
    printf("pitch=%.3f deg\n", tilt.pitch * 180 / pi);

    const float saved = tilt.pitch;
    ok = accel_tilt(0, 0, 0, &tilt);
    assert(!ok && tilt.pitch == saved);
    ok = accel_tilt(1, 0, 0, &tilt);
    assert(!ok); /* 俯仰 90 度附近的横滚奇异情况。 */
    ok = accel_tilt(NAN, 0, 1, &tilt);
    assert(!ok);
    puts("accel_tilt: PASS");
    return 0;
}

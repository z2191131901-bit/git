/*
 * 单轴互补滤波：用陀螺仪预测，用加速度计倾角校正。
 * 仅用于绕固定 X 轴的横滚教学模型，不能直接推广为任意三维运动。
 * gyro_rate 单位 rad/s；accel_angle、angle 单位 rad；dt、tau 单位 s。
 */
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

static const double PI = 3.14159265358979323846;

typedef struct
{
    float angle;
    int initialized;
} Complementary;

/* 把角度变成 [-pi,pi] 内的等效表示，避免跨越 180 度时走远路。 */
static double wrap_angle(double angle)
{
    return remainder(angle, 2.0 * PI);
}

/*
 * 使用 {0} 初始化。accel_valid 是上层对加速度倾角可信度的判断。
 * 还没初始化时必须有可信倾角；初始化后可以暂时仅用陀螺仪预测。
 * 不可信时允许 accel_angle 为 NaN，因为这次不会读取它参与计算。
 */
static int complementary_update(Complementary *state, float gyro_rate,
                                float accel_angle, int accel_valid,
                                float tau, float dt)
{
    if (state == NULL || !isfinite(gyro_rate) ||
        !isfinite(tau) || tau <= 0 || !isfinite(dt) || dt <= 0 ||
        (accel_valid && !isfinite(accel_angle)))
        return 0;

    if (!state->initialized)
    {
        if (!accel_valid)
            return 0;
        /* 第一条观测定义起始时刻，不再额外积分前一个未知时间段。 */
        state->angle = (float)wrap_angle(accel_angle);
        state->initialized = 1;
        return 1;
    }

    /* ① 预测：旧角度 + 角速度 × 时间。 */
    double next = wrap_angle((double)state->angle + (double)gyro_rate * dt);

    if (accel_valid)
    {
        /*
         * ② 求最短角度误差，再用 beta 比例修正。
         * beta = 1-alpha = dt/(tau+dt)，tau 越大，校正越慢。
         */
        const double beta = (double)dt / ((double)tau + dt);
        const double error = wrap_angle((double)accel_angle - next);
        next = wrap_angle(next + beta * error);
    }

    if (!isfinite(next))
        return 0;
    state->angle = (float)next;
    return 1;
}

int main(void)
{
    Complementary filter = {0};
    int ok = complementary_update(&filter, 0, NAN, 0, 0.5f, 0.01f);
    assert(!ok && !filter.initialized);

    ok = complementary_update(&filter, 0, 0, 1, 0.5f, 0.01f);
    assert(ok);
    ok = complementary_update(&filter, 0.1f, 0, 1, 0.5f, 0.01f);
    /* 预测 0.001 rad，再乘 alpha=0.5/0.51，得到约 0.000980392。 */
    assert(ok && fabsf(filter.angle - 0.000980392f) < 1e-7f);

    const float saved = filter.angle;
    ok = complementary_update(&filter, 0, 0, 1, 0.5f, 0);
    assert(!ok && filter.angle == saved);

    /* 暂无可靠加速度观测时，执行纯陀螺仪预测。 */
    ok = complementary_update(&filter, 0.1f, NAN, 0, 0.5f, 0.01f);
    assert(ok && fabsf(filter.angle - (saved + 0.001f)) < 1e-7f);

    /* +179 到 -179 度应只差 2 度，而不是差 -358 度。 */
    Complementary crossing = {(float)(179 * PI / 180), 1};
    ok = complementary_update(&crossing, 0, (float)(-179 * PI / 180),
                              1, 0.01f, 0.01f);
    assert(ok && fabs(fabs(crossing.angle) - PI) < 1e-5);
    printf("boundary result=%.3f deg (expected +/-180)\n",
           crossing.angle * 180 / PI);

    /* 陀螺仪零偏未去除时，互补滤波仍可能有稳态误差，不能声称完全消除。 */
    Complementary bias_case = {0, 1};
    for (int i = 0; i < 2000; ++i)
    {
        ok = complementary_update(&bias_case, 0.02f, 0, 1, 0.5f, 0.01f);
        assert(ok);
    }
    assert(fabsf(bias_case.angle - 0.01f) < 1e-5f);
    printf("bias steady error=%.6f rad\n", bias_case.angle);
    puts("complementary: PASS");
    return 0;
}

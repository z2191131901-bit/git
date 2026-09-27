/*
 * 幅值不变的 Clarke / Park 变换，角度单位 rad，正角度逆时针。
 * abc 相序为 a,b,c；保留零序，因而也能往返变换不平衡的三个瞬时值。
 * 此处只有代数运算，不负责电流采样、转子位置估计或电机闭环。
 */
#ifndef MOTOR_COORDINATES_H
#define MOTOR_COORDINATES_H
#include <math.h>
#include <stddef.h>

#define MOTOR_PI 3.14159265358979323846

typedef struct {
    double a, b, c; /* 三相瞬时量，必须采用相同单位，例如 A 或 V。 */
} MotorAbc;

typedef struct {
    double alpha, beta; /* 静止坐标系，alpha 与 a 轴重合。 */
    double zero;        /* 零序分量：(a+b+c)/3。 */
} MotorAb0;

typedef struct {
    double d, q; /* d 轴朝 theta，q 轴在 d 轴逆时针 90 度方向。 */
    double zero;
} MotorDq0;

/* 下列函数成功返回1、失败返回0；先在局部计算，失败不改调用者输出。 */
static inline int motor_clarke(MotorAbc x, MotorAb0 *out)
{
    MotorAb0 next;
    if (!out || !isfinite(x.a) || !isfinite(x.b) || !isfinite(x.c))
        return 0;
    next.alpha = (2.0 * x.a - x.b - x.c) / 3.0;
    next.beta = (x.b - x.c) / sqrt(3.0);
    next.zero = (x.a + x.b + x.c) / 3.0;
    if (!isfinite(next.alpha) || !isfinite(next.beta) || !isfinite(next.zero))
        return 0;
    *out = next;
    return 1;
}

static inline int motor_inverse_clarke(MotorAb0 x, MotorAbc *out)
{
    MotorAbc next;
    if (!out || !isfinite(x.alpha) || !isfinite(x.beta) || !isfinite(x.zero))
        return 0;
    next.a = x.alpha + x.zero;
    next.b = -0.5 * x.alpha + sqrt(3.0) * 0.5 * x.beta + x.zero;
    next.c = -0.5 * x.alpha - sqrt(3.0) * 0.5 * x.beta + x.zero;
    if (!isfinite(next.a) || !isfinite(next.b) || !isfinite(next.c))
        return 0;
    *out = next;
    return 1;
}

/* Park 是坐标轴旋转：矢量不变，改变它在旋转坐标系中的表示。 */
static inline int motor_park(MotorAb0 x, double theta, MotorDq0 *out)
{
    MotorDq0 next;
    double cs, sn;
    if (!out || !isfinite(theta) || !isfinite(x.alpha) ||
        !isfinite(x.beta) || !isfinite(x.zero))
        return 0;
    cs = cos(theta);
    sn = sin(theta);
    next.d = x.alpha * cs + x.beta * sn;
    next.q = -x.alpha * sn + x.beta * cs;
    next.zero = x.zero;
    if (!isfinite(next.d) || !isfinite(next.q))
        return 0;
    *out = next;
    return 1;
}

static inline int motor_inverse_park(MotorDq0 x, double theta, MotorAb0 *out)
{
    MotorAb0 next;
    double cs, sn;
    if (!out || !isfinite(theta) || !isfinite(x.d) ||
        !isfinite(x.q) || !isfinite(x.zero))
        return 0;
    cs = cos(theta);
    sn = sin(theta);
    next.alpha = x.d * cs - x.q * sn;
    next.beta = x.d * sn + x.q * cs;
    next.zero = x.zero;
    if (!isfinite(next.alpha) || !isfinite(next.beta))
        return 0;
    *out = next;
    return 1;
}
#endif

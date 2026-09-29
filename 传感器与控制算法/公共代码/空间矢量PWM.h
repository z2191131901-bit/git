/* 两电平SVPWM公共实现；第74篇给出手算和完整验证，第76篇FOC直接复用。 */
#ifndef MOTOR_SVPWM_H
#define MOTOR_SVPWM_H
#include <math.h>
#include "电机坐标变换.h"

typedef struct {
    double duty_a, duty_b, duty_c; /* 上桥臂导通占空比，范围[0,1]。 */
    double alpha, beta;           /* 限幅后真正要求调制器输出的电压[V]。 */
    double scale;                 /* 请求缩放比例，1为未削减，0~1为缩小。 */
    int limited;                  /* 请求是否超过本例采用的电压圆。 */
} PwmResult;

/* 只用于清除浮点计算在0或1附近的舍入误差，不代替矢量限幅。 */
static inline double duty_roundoff_clamp(double x)
{
    return fmax(0.0, fmin(1.0, x));
}

/* 成功返回1，非法输入返回0且不改*out。输出指针不能为NULL。 */
static inline int svpwm(double alpha, double beta, double bus, PwmResult *out)
{
    PwmResult next = {0};
    MotorAbc phase;
    double magnitude, maximum_voltage, common, highest, lowest;
    if (!out || !isfinite(alpha) || !isfinite(beta) ||
        !isfinite(bus) || bus < 1e-9)
        return 0;
    magnitude = hypot(alpha, beta);
    maximum_voltage = bus / sqrt(3.0);
    if (!isfinite(magnitude))
        return 0;

    /*
     * 按同一个比例缩小两个分量，保持电压矢量方向。
     * 不要分别截断alpha、beta，否则角度也会改变。
     */
    next.scale = magnitude > maximum_voltage ? maximum_voltage / magnitude : 1.0;
    next.limited = magnitude > maximum_voltage;
    next.alpha = alpha * next.scale;
    next.beta = beta * next.scale;
    if (!motor_inverse_clarke((MotorAb0){next.alpha, next.beta, 0.0}, &phase))
        return 0;

    /*
     * 三相加同一个common后，相间差值不变。
     * 把最大相与最小相的中点搬到0，可充分利用正负半母线范围。
     */
    highest = fmax(phase.a, fmax(phase.b, phase.c));
    lowest = fmin(phase.a, fmin(phase.b, phase.c));
    common = -0.5 * (highest + lowest);
    next.duty_a = 0.5 + (phase.a + common) / bus;
    next.duty_b = 0.5 + (phase.b + common) / bus;
    next.duty_c = 0.5 + (phase.c + common) / bus;

    /* 若超过允许舍入量，说明不能把它当成正常线性PWM结果。 */
    if (!isfinite(next.duty_a) || !isfinite(next.duty_b) || !isfinite(next.duty_c) ||
        next.duty_a < -1e-12 || next.duty_a > 1.0 + 1e-12 ||
        next.duty_b < -1e-12 || next.duty_b > 1.0 + 1e-12 ||
        next.duty_c < -1e-12 || next.duty_c > 1.0 + 1e-12)
        return 0;
    next.duty_a = duty_roundoff_clamp(next.duty_a);
    next.duty_b = duty_roundoff_clamp(next.duty_b);
    next.duty_c = duty_roundoff_clamp(next.duty_c);
    *out = next;
    return 1;
}

#endif

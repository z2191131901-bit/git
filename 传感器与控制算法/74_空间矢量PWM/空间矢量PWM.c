/*
 * 第74篇：两电平三相逆变器的空间矢量PWM。
 * 使用“逆Clarke + 最大最小值公共偏置”计算对称SVPWM的平均占空比。
 * 输入alpha/beta为相对负载中性点的电压矢量[V]，母线电压为Vdc[V]。
 * 本例保守地把矢量限制在线性内切圆Vdc/sqrt(3)中，不实现过调制。
 *
 * 这里只计算占空比，未实现定时器、死区、最小脉宽和电流环。
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "../公共代码/电机坐标变换.h"

typedef struct {
    double duty_a, duty_b, duty_c; /* 上桥臂导通占空比，范围[0,1]。 */
    double alpha, beta;           /* 限幅后真正要求调制器输出的电压[V]。 */
    double scale;                 /* 请求缩放比例，1为未削减，0~1为缩小。 */
    int limited;                  /* 请求是否超过本例采用的电压圆。 */
} PwmResult;

/* 只用于清除浮点计算在0或1附近的舍入误差，不代替矢量限幅。 */
static double duty_roundoff_clamp(double x)
{
    return fmax(0.0, fmin(1.0, x));
}

/* 成功返回1，非法输入返回0且不改*out。输出指针不能为NULL。 */
static int svpwm(double alpha, double beta, double bus, PwmResult *out)
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

int main(void)
{
    PwmResult result;
    double max_vector_error = 0.0;
    const double bus = 12.0;
    const double ratios[] = {0.0, 0.8, 1.0, 1.8};

    /*
     * 手算：alpha=6,beta=0 => 相电压(6,-3,-3)。
     * common=-1.5，故三相占空比为(0.875,0.125,0.125)。
     */
    assert(svpwm(6.0, 0.0, bus, &result));
    assert(fabs(result.duty_a - 0.875) < 1e-12);
    assert(fabs(result.duty_b - 0.125) < 1e-12 &&
           fabs(result.duty_c - 0.125) < 1e-12);
    assert(!result.limited);
    printf("SVPWM: hand duties=(%.3f, %.3f, %.3f)\n",
           result.duty_a, result.duty_b, result.duty_c);

    /*
     * 六个扇区都扫描，并从占空比反推实际平均电压。
     * 这里用“桥臂对母线中点的电压”，Clarke会消除公共偏置。
     */
    for (int level = 0; level < 4; ++level) {
        for (int degree = 0; degree < 360; ++degree) {
            double angle = degree * MOTOR_PI / 180.0;
            double request = ratios[level] * bus / sqrt(3.0);
            double expected = fmin(ratios[level], 1.0) * bus / sqrt(3.0);
            MotorAb0 reconstructed;
            assert(svpwm(request * cos(angle), request * sin(angle), bus, &result));
            MotorAbc pole_voltage = {
                bus * (result.duty_a - 0.5),
                bus * (result.duty_b - 0.5),
                bus * (result.duty_c - 0.5)
            };
            assert(motor_clarke(pole_voltage, &reconstructed));
            double error = hypot(reconstructed.alpha - expected * cos(angle),
                                 reconstructed.beta - expected * sin(angle));
            max_vector_error = fmax(max_vector_error, error);
            assert(result.duty_a >= 0.0 && result.duty_a <= 1.0);
            assert(result.duty_b >= 0.0 && result.duty_b <= 1.0);
            assert(result.duty_c >= 0.0 && result.duty_c <= 1.0);
            assert(fabs((result.duty_a - result.duty_b) * bus -
                        (1.5 * result.alpha - sqrt(3.0) * 0.5 * result.beta)) < 1e-11);
            if (level == 3)
                assert(result.limited && fabs(result.scale - 1.0 / 1.8) < 1e-12);
        }
    }
    assert(max_vector_error < 1e-11);

    /* 30度处落在内切圆与六边形的切点，占空比应为(1,0.5,0)。 */
    assert(svpwm(bus / 2.0, bus / (2.0 * sqrt(3.0)), bus, &result));
    assert(fabs(result.duty_a - 1.0) < 1e-12);
    assert(fabs(result.duty_b - 0.5) < 1e-12 && fabs(result.duty_c) < 1e-12);

    /* 零请求应给三路相同的50%占空比，相间平均电压为0。 */
    assert(svpwm(0.0, 0.0, bus, &result));
    assert(result.duty_a == 0.5 && result.duty_b == 0.5 && result.duty_c == 0.5);
    assert(!svpwm(0.0, 0.0, 0.0, &result));
    assert(!svpwm(NAN, 0.0, bus, &result));
    assert(result.duty_a == 0.5);
    assert(!svpwm(0.0, 0.0, bus, NULL));

    printf("SVPWM: reconstructed_vector_max_error=%.3g V; circle_limit=%.9f V\n",
           max_vector_error, bus / sqrt(3.0));
    puts("PASS: 1440 vectors, sector boundaries, radial limiting and invalid inputs.");
    return 0;
}

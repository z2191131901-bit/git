/*
 * 第73篇：Clarke / Park 与逆变换。
 * 算法函数在公共中文头文件中，第74、75篇复用同一坐标约定。
 * 本文件用手算、三相旋转信号和不平衡信号检查公式。
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "../公共代码/电机坐标变换.h"

int main(void)
{
    MotorAb0 stationary, restored;
    MotorDq0 rotating;
    MotorAbc phases;
    double max_roundtrip = 0.0, max_power_error = 0.0;
    const MotorAbc hand = {2.0, -1.0, -1.0};

    /* 手算：abc=(2,-1,-1) => alpha=2,beta=0,zero=0。 */
    assert(motor_clarke(hand, &stationary));
    assert(fabs(stationary.alpha - 2.0) < 1e-12);
    assert(fabs(stationary.beta) < 1e-12 && fabs(stationary.zero) < 1e-12);

    /* 把坐标轴转90度，同一个矢量在新坐标中变成 d=0,q=-2。 */
    assert(motor_park(stationary, MOTOR_PI / 2.0, &rotating));
    assert(fabs(rotating.d) < 1e-12 && fabs(rotating.q + 2.0) < 1e-12);

    for (int k = 0; k < 721; ++k) {
        double theta = -MOTOR_PI + k * MOTOR_PI / 180.0;
        double current_peak = 3.0;
        MotorAbc balanced = {
            current_peak * cos(theta),
            current_peak * cos(theta - 2.0 * MOTOR_PI / 3.0),
            current_peak * cos(theta + 2.0 * MOTOR_PI / 3.0)
        };

        /*
         * 正确相序、正确角度下，旋转三相电流变成恒定d轴量。
         * 这里的3是相电流峰值，不是有效值。
         */
        assert(motor_clarke(balanced, &stationary));
        assert(motor_park(stationary, theta, &rotating));
        assert(fabs(rotating.d - current_peak) < 1e-12);
        assert(fabs(rotating.q) < 1e-12);

        /* 角度估计多出0.1 rad时，q应为 -3*sin(0.1)，核对误差符号。 */
        assert(motor_park(stationary, theta + 0.1, &rotating));
        assert(fabs(rotating.q + current_peak * sin(0.1)) < 1e-12);

        /*
         * 非平衡输入也要能恢复，因此不能丢掉zero。
         * 相量平方和在幅值不变变换下有3/2系数，不是数值直接相等。
         */
        MotorAbc unbalanced = {1.0 + cos(theta), -0.4, 0.7 * sin(2.0 * theta)};
        assert(motor_clarke(unbalanced, &stationary));
        assert(motor_park(stationary, theta, &rotating));
        assert(motor_inverse_park(rotating, theta, &restored));
        assert(motor_inverse_clarke(restored, &phases));
        double error = fmax(fabs(phases.a - unbalanced.a),
                       fmax(fabs(phases.b - unbalanced.b),
                            fabs(phases.c - unbalanced.c)));
        double abc_square = unbalanced.a * unbalanced.a +
                            unbalanced.b * unbalanced.b + unbalanced.c * unbalanced.c;
        double dq_square = 1.5 * (rotating.d * rotating.d + rotating.q * rotating.q) +
                           3.0 * rotating.zero * rotating.zero;
        max_roundtrip = fmax(max_roundtrip, error);
        max_power_error = fmax(max_power_error, fabs(abc_square - dq_square));
    }
    assert(max_roundtrip < 1e-12 && max_power_error < 1e-11);

    /* 三相同为1是纯零序：alpha=beta=0，zero=1。 */
    assert(motor_clarke((MotorAbc){1.0, 1.0, 1.0}, &stationary));
    assert(stationary.alpha == 0.0 && stationary.beta == 0.0 && stationary.zero == 1.0);

    /* 非法角度或输入不能写坏上一份结果。 */
    MotorDq0 saved = rotating;
    assert(!motor_park(stationary, NAN, &rotating));
    assert(rotating.d == saved.d && rotating.q == saved.q);
    assert(!motor_clarke((MotorAbc){NAN, 0.0, 0.0}, &stationary));
    assert(stationary.zero == 1.0);
    assert(!motor_inverse_clarke(stationary, NULL));

    printf("coordinates: roundtrip_max=%.3g square_identity_max=%.3g\n",
           max_roundtrip, max_power_error);
    puts("PASS: hand values, phase sequence, zero sequence and inverse transforms.");
    return 0;
}

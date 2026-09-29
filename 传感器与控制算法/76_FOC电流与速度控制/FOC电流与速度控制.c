/*
 * 第76篇：带理想角度反馈的PMSM FOC电流环与速度环。
 * 控制链：三相电流 -> Clarke/Park -> dq电流PI与解耦 ->
 *         逆Park -> SVPWM -> 平均电压 -> 电机 -> 下一次测量。
 *
 * 电机在静止alpha/beta坐标下积分，控制器在dq坐标下计算。
 * 这样仿真不仅验证PI，也检查坐标方向和PWM链是否接对。
 * 本例是平均模型教学，没有开关纹波、编码器误差或真实驱动硬件。
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "../公共代码/空间矢量PWM.h"
#include "../公共代码/PID核心.h"

/* 表贴式PMSM参数：Ld=Lq，电流使用幅值不变变换后的峰值约定。 */
static const double RESISTANCE = 0.4;  /* 相电阻[ohm]。 */
static const double INDUCTANCE = 0.001; /* d/q轴电感[H]。 */
static const double FLUX = 0.02;       /* 永磁体磁链[Wb]。 */
static const double POLE_PAIRS = 4.0;  /* 极对数；不是磁极总数。 */
static const double INERTIA = 0.0005;  /* 转动惯量[kg*m^2]。 */
static const double FRICTION = 0.0001; /* 黏性摩擦[N*m/(rad/s)]。 */
static const double DT = 0.0001;       /* 电流环周期[s]，10kHz。 */

typedef struct {
    double alpha, beta; /* 静止坐标电流[A]。 */
    double speed;       /* 机械角速度[rad/s]。 */
    double theta;       /* 电角度[rad]，运行状态保持在[0,2*pi)。 */
} MotorState;

typedef struct {
    double integral_d, integral_q; /* 已乘Ki的电流积分贡献[V]。 */
} CurrentLoop;

/* 电压与转矩作用于真实对象；锁轴模式只积分电流，外部夹具承担转矩。 */
static MotorState derivative(MotorState x, double va, double vb,
                             double load, int locked)
{
    const double electrical_speed = POLE_PAIRS * x.speed;
    const double iq = -x.alpha * sin(x.theta) + x.beta * cos(x.theta);
    const double torque = 1.5 * POLE_PAIRS * FLUX * iq;
    MotorState dx;
    dx.alpha = (va - RESISTANCE * x.alpha +
                electrical_speed * FLUX * sin(x.theta)) / INDUCTANCE;
    dx.beta = (vb - RESISTANCE * x.beta -
               electrical_speed * FLUX * cos(x.theta)) / INDUCTANCE;
    dx.speed = locked ? 0.0 : (torque - load - FRICTION * x.speed) / INERTIA;
    dx.theta = locked ? 0.0 : electrical_speed;
    return dx;
}

static MotorState add_scaled(MotorState x, MotorState dx, double h)
{
    MotorState result = {
        x.alpha + h * dx.alpha, x.beta + h * dx.beta,
        x.speed + h * dx.speed, x.theta + h * dx.theta
    };
    return result;
}

/*
 * 一个控制周期内alpha/beta平均电压保持不变。
 * 用RK4子步积分连续电机方程；子步数只改变仿真精度，不改变控制周期。
 */
static void motor_step(MotorState *x, double va, double vb, double load,
                       int locked, int substeps)
{
    const double h = DT / substeps;
    for (int j = 0; j < substeps; ++j) {
        MotorState a = derivative(*x, va, vb, load, locked);
        MotorState b = derivative(add_scaled(*x, a, h / 2.0), va, vb, load, locked);
        MotorState c = derivative(add_scaled(*x, b, h / 2.0), va, vb, load, locked);
        MotorState d = derivative(add_scaled(*x, c, h), va, vb, load, locked);
        x->alpha += h / 6.0 * (a.alpha + 2.0*b.alpha + 2.0*c.alpha + d.alpha);
        x->beta += h / 6.0 * (a.beta + 2.0*b.beta + 2.0*c.beta + d.beta);
        x->speed += h / 6.0 * (a.speed + 2.0*b.speed + 2.0*c.speed + d.speed);
        x->theta += h / 6.0 * (a.theta + 2.0*b.theta + 2.0*c.theta + d.theta);
    }
    x->theta = fmod(x->theta, 2.0 * MOTOR_PI);
    if (x->theta < 0.0)
        x->theta += 2.0 * MOTOR_PI;
    assert(isfinite(x->alpha) && isfinite(x->beta) && isfinite(x->speed));
}

/* 从占空比重建实际平均桥臂电压，而不是把未限幅的PI请求直接给对象。 */
static MotorAb0 actual_voltage(PwmResult pwm, double bus)
{
    MotorAbc pole = {
        (pwm.duty_a - 0.5) * bus,
        (pwm.duty_b - 0.5) * bus,
        (pwm.duty_c - 0.5) * bus
    };
    MotorAb0 voltage;
    assert(motor_clarke(pole, &voltage));
    return voltage;
}

/*
 * 一次电流环更新；theta必须是当前电角度，speed是机械角速度。
 * 参考与测量单位A，bus为V；积分初始设0。
 * antiwindup非零时启用矢量限幅后的反算；为0仅用于对照实验。
 * 成功返回1；非法输入返回0且不改积分或PWM结果。
 * 输出不得与state重叠。调度周期在本例固定为DT。
 */
static int foc_step(CurrentLoop *state, MotorAbc measured, double theta,
                    double speed, double id_reference, double iq_reference,
                    double bus, int antiwindup, PwmResult *output)
{
    const double bandwidth = 1800.0; /* 电流环设计角频率参数[rad/s]。 */
    const double kp = INDUCTANCE * bandwidth; /* 1.8 V/A。 */
    const double ki = RESISTANCE * bandwidth; /* 720 V/(A*s)。 */
    const double kaw = antiwindup ? 500.0 : 0.0; /* 反算系数[1/s]。 */
    MotorAb0 stationary, voltage_request;
    MotorDq0 current, applied;
    PwmResult pwm;
    CurrentLoop next;
    if (!state || !output || !isfinite(theta) || !isfinite(speed) ||
        !isfinite(id_reference) || !isfinite(iq_reference) ||
        !isfinite(state->integral_d) || !isfinite(state->integral_q))
        return 0;
    if (!motor_clarke(measured, &stationary) ||
        !motor_park(stationary, theta, &current))
        return 0;

    const double we = POLE_PAIRS * speed;
    const double ed = id_reference - current.d;
    const double eq = iq_reference - current.q;
    /*
     * 对象方程：
     * vd=R*id+L*id_dot-we*L*iq
     * vq=R*iq+L*iq_dot+we*(L*id+flux)
     * 下面的前馈抵消转速导致的交叉耦合与反电动势。
     */
    MotorDq0 raw = {
        kp * ed + state->integral_d - we * INDUCTANCE * current.q,
        kp * eq + state->integral_q + we * (INDUCTANCE * current.d + FLUX),
        0.0
    };
    if (!motor_inverse_park(raw, theta, &voltage_request) ||
        !svpwm(voltage_request.alpha, voltage_request.beta, bus, &pwm))
        return 0;

    MotorAb0 realized = actual_voltage(pwm, bus);
    if (!motor_park(realized, theta, &applied))
        return 0;
    /*
     * PWM的圆形限幅同时作用于两个轴。
     * 用实际电压减原始请求反算，才能包含前馈占用的母线余量。
     * 积分在本拍输出计算之后更新，供下一拍使用。
     */
    next.integral_d = state->integral_d + DT * (ki * ed + kaw * (applied.d - raw.d));
    next.integral_q = state->integral_q + DT * (ki * eq + kaw * (applied.q - raw.q));
    if (!isfinite(next.integral_d) || !isfinite(next.integral_q))
        return 0;
    *state = next;
    *output = pwm;
    return 1;
}

/* 模拟理想三相电流传感器；输出后由控制器重新做Clarke/Park。 */
static MotorAbc measured_currents(MotorState x)
{
    MotorAbc phases;
    assert(motor_inverse_clarke((MotorAb0){x.alpha, x.beta, 0.0}, &phases));
    return phases;
}

static void check_model_and_current(void)
{
    MotorState x = {0};
    CurrentLoop loop = {0};
    PwmResult pwm = {0};

    /* 锁轴R-L对象的独立解析解，用于检查连续模型积分。 */
    for (int k = 0; k < 100; ++k)
        motor_step(&x, 0.0, 1.0, 0.0, 1, 4);
    double expected = (1.0 / RESISTANCE) * (1.0 - exp(-RESISTANCE * 100 * DT / INDUCTANCE));
    assert(fabs(x.beta - expected) < 1e-9 && x.alpha == 0.0);

    /* 零转速、零初始电流，iq目标2A：第一拍vq=1.8*2=3.6V。 */
    x = (MotorState){0};
    assert(foc_step(&loop, measured_currents(x), 0.0, 0.0, 0.0, 2.0, 24.0, 1, &pwm));
    MotorAb0 voltage = actual_voltage(pwm, 24.0);
    assert(fabs(voltage.alpha) < 1e-12 && fabs(voltage.beta - 3.6) < 1e-12);
    assert(fabs(loop.integral_q - 0.144) < 1e-12);

    /* 无效母线电压不允许污染已经建立的积分和输出。 */
    double saved_i = loop.integral_q, saved_duty = pwm.duty_a;
    assert(!foc_step(&loop, measured_currents(x), 0.0, 0.0, 0.0, 2.0, 0.0, 1, &pwm));
    assert(loop.integral_q == saved_i && pwm.duty_a == saved_duty);

    loop = (CurrentLoop){0};
    for (int k = 0; k < 2000; ++k) {
        double reference = k < 1000 ? 2.0 : -1.0;
        assert(foc_step(&loop, measured_currents(x), x.theta, 0.0, 0.0,
                        reference, 24.0, 1, &pwm));
        voltage = actual_voltage(pwm, 24.0);
        motor_step(&x, voltage.alpha, voltage.beta, 0.0, 1, 4);
        if (k == 999)
            assert(fabs(x.beta - 2.0) < 1e-6);
    }
    assert(fabs(x.beta + 1.0) < 1e-6 && fabs(x.alpha) < 1e-10);
    printf("FOC: locked rotor final id=%.9f A iq=%.9f A\n", x.alpha, x.beta);
}

typedef struct {
    double recovery_iae; /* 母线恢复后电流误差的离散积分[A*s]。 */
    double peak;         /* 母线恢复后的最大绝对q轴电流[A]。 */
    int limited;
} Recovery;

static Recovery saturation_test(int antiwindup)
{
    MotorState x = {0};
    CurrentLoop loop = {0};
    Recovery result = {0};
    for (int k = 0; k < 3000; ++k) {
        double bus = k < 1000 ? 1.0 : 24.0;
        double reference = k < 1000 ? 4.0 : 1.0;
        PwmResult pwm;
        assert(foc_step(&loop, measured_currents(x), 0.0, 0.0, 0.0, reference,
                        bus, antiwindup, &pwm));
        result.limited += pwm.limited;
        MotorAb0 voltage = actual_voltage(pwm, bus);
        motor_step(&x, voltage.alpha, voltage.beta, 0.0, 1, 4);
        if (k >= 1000) {
            result.recovery_iae += fabs(reference - x.beta) * DT;
            result.peak = fmax(result.peak, fabs(x.beta));
        }
    }
    assert(fabs(x.beta - 1.0) < 1e-5);
    return result;
}

typedef struct {
    double final_speed, final_iq, final_id;
    double maximum_current, maximum_iq_reference;
    double before_load, after_load, recovery_iae;
    int voltage_limited;
} SpeedResult;

static SpeedResult speed_test(int substeps)
{
    MotorState x = {0};
    CurrentLoop current_loop = {0};
    PidState speed_loop = {0};
    /* 速度PI输出是iq参考[A]，限于±4A；每10个电流周期更新一次。 */
    /* 理想电流跟踪下，速度对象为Kt/(J*s+B)。按二阶目标配速度PI。 */
    const double kt = 1.5 * POLE_PAIRS * FLUX;
    const double wn = 20.0, zeta = 0.9; /* 自然角频率[rad/s]与阻尼比。 */
    const PidConfig cfg = {
        (2.0 * zeta * wn * INERTIA - FRICTION) / kt,
        INERTIA * wn * wn / kt, 0.0, 0.0, -4.0, 4.0
    };
    double iq_reference = 0.0;
    SpeedResult result = {0};
    for (int k = 0; k < 25000; ++k) {
        double reference = k < 14000 ? 80.0 : -60.0; /* 机械rad/s。 */
        double load = k < 6000 ? 0.0 : k < 14000 ? 0.06 : -0.04; /* 有符号N*m。 */
        if (k % 10 == 0)
            assert(pid_update(&speed_loop, &cfg, reference, x.speed,
                              0.0, 10.0 * DT, &iq_reference));
        assert(fabs(iq_reference) <= 4.0);
        result.maximum_iq_reference = fmax(result.maximum_iq_reference, fabs(iq_reference));

        PwmResult pwm;
        assert(foc_step(&current_loop, measured_currents(x), x.theta, x.speed,
                        0.0, iq_reference, 24.0, 1, &pwm));
        result.voltage_limited += pwm.limited;
        MotorAb0 voltage = actual_voltage(pwm, 24.0);
        motor_step(&x, voltage.alpha, voltage.beta, load, 0, substeps);
        result.maximum_current = fmax(result.maximum_current, hypot(x.alpha, x.beta));
        if (k == 5999)
            result.before_load = x.speed;
        if (k == 13999)
            result.after_load = x.speed;
        if (k >= 6000 && k < 14000)
            result.recovery_iae += fabs(80.0 - x.speed) * DT;
    }
    MotorDq0 current;
    assert(motor_park((MotorAb0){x.alpha, x.beta, 0.0}, x.theta, &current));
    result.final_speed = x.speed;
    result.final_iq = current.q;
    result.final_id = current.d;
    return result;
}

int main(void)
{
    check_model_and_current();
    Recovery without = saturation_test(0);
    Recovery with = saturation_test(1);
    assert(with.limited > 0 && without.limited > 0);
    assert(with.recovery_iae < 0.2 * without.recovery_iae);
    assert(with.peak < without.peak);
    printf("FOC: voltage recovery IAE without=%.9f with=%.9f A*s; peak=%.6f/%.6f A\n",
           without.recovery_iae, with.recovery_iae, without.peak, with.peak);

    SpeedResult result = speed_test(4);
    SpeedResult refined = speed_test(8);
    printf("FOC: RK4 refinement speed_difference=%.9g rad/s IAE_difference=%.9g rad\n",
           result.final_speed - refined.final_speed,
           result.recovery_iae - refined.recovery_iae);
    /* 子步加密是检查连续对象积分误差，不是改变控制器采样周期。 */
    assert(fabs(result.final_speed - refined.final_speed) < 1e-6);
    assert(fabs(result.recovery_iae - refined.recovery_iae) < 1e-5);
    assert(fabs(result.before_load - 80.0) < 0.1);
    assert(fabs(result.after_load - 80.0) < 0.1);
    assert(fabs(result.final_speed + 60.0) < 0.1);
    assert(fabs(result.final_id) < 0.01);

    /*
     * 稳态机械力矩平衡：Kt*iq=load+B*speed。
     * 最后一段load=-0.04，speed=-60，所以iq理论值约-0.383333A。
     */
    double expected_iq = (-0.04 + FRICTION * -60.0) / (1.5 * POLE_PAIRS * FLUX);
    assert(fabs(result.final_iq - expected_iq) < 0.005);
    printf("FOC: speed before_load=%.6f after_load=%.6f final=%.6f rad/s\n",
           result.before_load, result.after_load, result.final_speed);
    printf("FOC: final id=%.9f iq=%.9f A; peak_current=%.6f A iq_ref_peak=%.6f A limited=%d\n",
           result.final_id, result.final_iq, result.maximum_current,
           result.maximum_iq_reference, result.voltage_limited);
    puts("PASS: current loop, vector anti-windup, speed/load/reversal and RK4 refinement.");
    return 0;
}

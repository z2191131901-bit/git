/*
 * 第75篇：同步旋转坐标系数字锁相环（SRF-PLL）。
 * 输入alpha/beta是正交分量，可由平衡三相信号经Clarke得到。
 * 不是直接接一条单相正弦的PLL；单相应用需先产生正交信号。
 *
 * 鉴相：q/幅值 = sin(真实相位-估计相位)。
 * 环路：PI产生角频率修正，再积分角频率生成下一拍相位。
 * 状态theta在调用前对应k时刻，调用后对应k+1时刻。
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "../公共代码/电机坐标变换.h"

typedef struct {
    double dt;             /* 采样周期[s]。 */
    double kp;             /* 归一化相位误差到角频率修正的比例系数[1/s]。 */
    double ki;             /* 积分系数[1/s^2]，rad按无量纲处理。 */
    double nominal_omega;  /* 名义角频率[rad/s]。 */
    double minimum_omega;  /* 估计角频率下限[rad/s]，本例只跟踪正频率。 */
    double maximum_omega;  /* 估计角频率上限[rad/s]。 */
    double minimum_amplitude; /* 输入幅值门限，与alpha/beta同单位。 */
} PllConfig;

typedef struct {
    PllConfig cfg;
    double theta;       /* [0,2*pi)的相位[rad]，每周期更新一次。 */
    double integral;    /* 已乘Ki的积分贡献[rad/s]。 */
    double omega;       /* 本拍用于相位积分的角频率[rad/s]。 */
    double detector;    /* 本拍归一化q轴误差；弱信号时设0，仅表示未鉴相。 */
    int ready;
} Pll;

static double wrap_angle(double angle)
{
    double wrapped = fmod(angle, 2.0 * MOTOR_PI);
    return wrapped < 0.0 ? wrapped + 2.0 * MOTOR_PI : wrapped;
}

/* 配置和初始角度合法时清空积分；失败保留原结构。 */
static int pll_init(Pll *p, PllConfig cfg, double initial_angle)
{
    Pll next = {0};
    if (!p || !isfinite(initial_angle) || !isfinite(cfg.dt) || cfg.dt <= 0.0 ||
        !isfinite(cfg.kp) || cfg.kp < 0.0 || !isfinite(cfg.ki) || cfg.ki < 0.0 ||
        !isfinite(cfg.nominal_omega) || !isfinite(cfg.minimum_omega) ||
        !isfinite(cfg.maximum_omega) || cfg.minimum_omega <= 0.0 ||
        cfg.minimum_omega >= cfg.maximum_omega ||
        cfg.nominal_omega < cfg.minimum_omega || cfg.nominal_omega > cfg.maximum_omega ||
        !isfinite(cfg.minimum_amplitude) || cfg.minimum_amplitude <= 0.0)
        return 0;
    /* 范围检查只避免明显的采样不足，不能代替完整离散稳定性设计。 */
    if (!isfinite(cfg.maximum_omega * cfg.dt) ||
        cfg.maximum_omega * cfg.dt >= 0.5 * MOTOR_PI)
        return 0;
    next.cfg = cfg;
    next.theta = wrap_angle(initial_angle);
    next.omega = cfg.nominal_omega;
    next.ready = 1;
    *p = next;
    return 1;
}

/*
 * 返回1：信号幅值足够，完成一次鉴相；不表示已经锁定。
 * 返回0：幅值不足，冻结积分和频率，按最近频率继续走相位。
 * 返回-1：参数/计算无效，整个状态不变。
 * 弱信号期间继续走相位是保持运行策略，不等于知道真实相位。
 */
static int pll_step(Pll *p, double alpha, double beta)
{
    Pll next;
    double amplitude;
    int measured;
    if (!p || !p->ready || !isfinite(alpha) || !isfinite(beta))
        return -1;
    amplitude = hypot(alpha, beta);
    if (!isfinite(amplitude))
        return -1;
    next = *p;
    measured = amplitude >= p->cfg.minimum_amplitude;
    next.detector = 0.0;

    if (measured) {
        MotorDq0 dq;
        /* 先归一化再Park，使PI参数不随输入幅值同比变化。 */
        MotorAb0 normalized = {alpha / amplitude, beta / amplitude, 0.0};
        if (!motor_park(normalized, p->theta, &dq))
            return -1;
        next.detector = dq.q;
        double increment = p->cfg.ki * next.detector * p->cfg.dt;
        double candidate = p->integral + increment;
        double raw = p->cfg.nominal_omega + p->cfg.kp * next.detector + candidate;
        if (!isfinite(candidate) || !isfinite(raw))
            return -1;

        /*
         * 频率达到上限却还想向上积分，或达到下限还想向下积分：
         * 拒绝本次积分，避免长时间积累后难以恢复。
         * 反方向积分仍允许。这里是条件积分抗饱和。
         */
        if ((raw > p->cfg.maximum_omega && increment > 0.0) ||
            (raw < p->cfg.minimum_omega && increment < 0.0)) {
            candidate = p->integral;
            raw = p->cfg.nominal_omega + p->cfg.kp * next.detector + candidate;
        }
        if (!isfinite(raw))
            return -1;
        next.integral = candidate;
        next.omega = fmax(p->cfg.minimum_omega, fmin(p->cfg.maximum_omega, raw));
    }

    /* 必须使用当前误差算出的频率推进下一拍，调用者不能再重复积分。 */
    double advanced = p->theta + next.omega * p->cfg.dt;
    if (!isfinite(advanced))
        return -1;
    next.theta = wrap_angle(advanced);
    *p = next;
    return measured;
}

static PllConfig example_config(void)
{
    /*
     * 小误差线性化特征式s^2+Kp*s+Ki。
     * 选择自然角频率wn=2*pi*10、阻尼比sqrt(1/2)。
     * 这里的10Hz不是直接声明闭环-3dB带宽。
     */
    const double wn = 2.0 * MOTOR_PI * 10.0;
    PllConfig cfg = {
        0.0005, sqrt(2.0) * wn, wn * wn,
        2.0 * MOTOR_PI * 50.0,
        2.0 * MOTOR_PI * 40.0, 2.0 * MOTOR_PI * 70.0, 0.05
    };
    return cfg;
}

static void check_boundaries(void)
{
    Pll p;
    PllConfig cfg = example_config();
    assert(pll_init(&p, cfg, 0.0));

    /* 真实输入领先估计0.1rad，鉴相误差应为sin(0.1)，推动频率增加。 */
    double error = sin(0.1);
    assert(pll_step(&p, cos(0.1), sin(0.1)) == 1);
    double expected = cfg.nominal_omega + cfg.kp * error + cfg.ki * error * cfg.dt;
    assert(fabs(p.omega - expected) < 1e-12);
    assert(fabs(p.theta - expected * cfg.dt) < 1e-12);

    double old_theta = p.theta, old_omega = p.omega, old_integral = p.integral;
    assert(pll_step(&p, 0.0, 0.0) == 0);
    assert(p.omega == old_omega && p.integral == old_integral);
    assert(fabs(p.theta - wrap_angle(old_theta + old_omega * cfg.dt)) < 1e-12);
    old_theta = p.theta;
    assert(pll_step(&p, NAN, 0.0) == -1 && p.theta == old_theta);

    /* 单独制造持续饱和，验证积分没有累积；这些输入不代表真实波形。 */
    cfg.kp = 400.0;
    cfg.maximum_omega = 2.0 * MOTOR_PI * 60.0;
    assert(pll_init(&p, cfg, 0.0));
    for (int k = 0; k < 20; ++k) {
        assert(pll_step(&p, -sin(p.theta), cos(p.theta)) == 1); /* q=+1。 */
        assert(p.omega == cfg.maximum_omega && p.integral == 0.0);
    }
    assert(pll_step(&p, cos(p.theta), sin(p.theta)) == 1); /* q=0，退出饱和。 */
    assert(fabs(p.omega - cfg.nominal_omega) < 1e-12);
    assert(pll_step(&p, sin(p.theta), -cos(p.theta)) == 1); /* q=-1。 */
    assert(p.omega == cfg.minimum_omega && fabs(p.integral) < 1e-12);
    cfg.dt = 0.0;
    old_theta = p.theta;
    assert(!pll_init(&p, cfg, 0.0) && p.theta == old_theta);
}

int main(void)
{
    Pll p;
    PllConfig cfg = example_config();
    double actual_theta = 0.7;
    double maximum_phase[5] = {0}, maximum_hz[5] = {0};
    int dropout_count = 0;
    check_boundaries();
    assert(pll_init(&p, cfg, 0.0));

    /*
     * 0s：初始相位差0.7rad。
     * 0.5s：真实相位再跳变0.4rad。
     * 1s：频率50变55Hz。1.5s：幅值2变0.4。
     * 2~2.1s：输入丢失，但真实频率仍为55Hz。
     */
    for (int k = 0; k < 5000; ++k) {
        double frequency = k < 2000 ? 50.0 : 55.0;
        double amplitude = k < 3000 ? 2.0 : 0.4;
        if (k == 1000)
            actual_theta = wrap_angle(actual_theta + 0.4);
        if (k >= 4000 && k < 4200)
            amplitude = 0.0;

        /* 比较同一时刻k的相位，不能拿真实k与估计k+1错位比较。 */
        double phase_error = atan2(sin(actual_theta - p.theta), cos(actual_theta - p.theta));
        MotorAbc phases = {
            amplitude * cos(actual_theta),
            amplitude * cos(actual_theta - 2.0 * MOTOR_PI / 3.0),
            amplitude * cos(actual_theta + 2.0 * MOTOR_PI / 3.0)
        };
        MotorAb0 stationary;
        assert(motor_clarke(phases, &stationary));
        int status = pll_step(&p, stationary.alpha, stationary.beta);
        assert(status == (amplitude == 0.0 ? 0 : 1));
        dropout_count += status == 0;
        assert(p.theta >= 0.0 && p.theta < 2.0 * MOTOR_PI);

        /* 每次变化后留出恢复时间，在明确窗口内统计最大误差。 */
        int window = k >= 800 && k < 1000 ? 0 :
                     k >= 1800 && k < 2000 ? 1 :
                     k >= 2800 && k < 3000 ? 2 :
                     k >= 3800 && k < 4000 ? 3 :
                     k >= 4600 && k < 5000 ? 4 : -1;
        if (window >= 0) {
            maximum_phase[window] = fmax(maximum_phase[window], fabs(phase_error));
            maximum_hz[window] = fmax(maximum_hz[window],
                                     fabs(p.omega / (2.0 * MOTOR_PI) - frequency));
        }
        actual_theta = wrap_angle(actual_theta + 2.0 * MOTOR_PI * frequency * cfg.dt);
    }
    assert(dropout_count == 200);
    for (int j = 0; j < 5; ++j) {
        assert(maximum_phase[j] < 0.003 && maximum_hz[j] < 0.02);
        printf("PLL: window%d max_phase=%.9g rad max_frequency=%.9g Hz\n",
               j, maximum_phase[j], maximum_hz[j]);
    }
    printf("PLL: final_frequency=%.9f Hz holdover_samples=%d\n",
           p.omega / (2.0 * MOTOR_PI), dropout_count);
    puts("PASS: phase/frequency steps, amplitude change, dropout and anti-windup.");
    return 0;
}

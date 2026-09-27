/*
 * 第70篇：比例谐振 PR 控制（有限带宽的准 PR）
 * 编译与运行：在仓库根目录执行讲解.md中的验证命令。
 *
 * 对象：tau * y_dot + y = u，y、u 都使用归一化量。
 * 控制器：C(s) = kp + 2*kr*wc*s / (s*s + 2*wc*s + w0*w0)。
 * kr 是中心频率处谐振支路的增益；它不是所有资料中的 Ki。
 * 采用中心频率预畸变的 Tustin 离散，避免谐振中心发生映射偏移。
 *
 * 这是没有执行器饱和的线性教学模型，没有实现输出限幅或抗饱和。
 * 主函数中的幅度 assert 只是检查本例数值，不会修改控制输出。
 */
#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>

#define PI 3.14159265358979323846

typedef struct {
    double kp;               /* 比例增益，单位为“输出单位/误差单位”。 */
    double b0, a1, a2;        /* 归一化离散系数；谐振分子为 b0*(1-z^-2)。 */
    double e1, e2;           /* 上一次、上上次误差。 */
    double v1, v2;           /* 上一次、上上次谐振支路输出。 */
    int ready;               /* 配置成功以后为 1。 */
} QuasiPr;

typedef struct {
    double rms_error;        /* 最后两秒误差的均方根。 */
    double predicted_rms;    /* 离散闭环频率响应独立预测的均方根。 */
    double peak_output;      /* 包含启动过程的控制输出绝对峰值。 */
} Result;

/*
 * 参数：w0 是中心角频率 [rad/s]，wc 是阻尼带宽参数 [rad/s]，
 * dt 是周期 [s]；kp、kr 非负，kr 为 0 时退化为 P 控制器。
 * 配置成功时清空历史量，失败时保留原控制器。
 * 不要在运行中每次调用本函数，否则谐振记忆会被不断清零。
 */
static int pr_configure(QuasiPr *controller, double kp, double kr,
                        double w0, double wc, double dt)
{
    QuasiPr next = {0};
    double k, d0;
    if (controller == NULL || !isfinite(kp) || !isfinite(kr) ||
        !isfinite(w0) || !isfinite(wc) || !isfinite(dt) ||
        kp < 0.0 || kr < 0.0 || w0 <= 0.0 || wc <= 0.0 || dt <= 0.0)
        return 0;

    /*
     * 本例主动把配置范围限制到中心频率低于采样频率的 40%。
     * 这只是数值/采样范围检查，不等于保证闭环稳定或工程采样充足。
     * 仿真实际为 50 Hz / 5000 Hz = 1%，每周期有 100 个样本。
     */
    if (!isfinite(w0 * dt) || w0 * dt >= 0.8 * PI)
        return 0;

    /* 普通 Tustin 使用 k=2/dt；预畸变使指定 w0 精确映射到自己。 */
    k = w0 / tan(w0 * dt / 2.0);
    d0 = k * k + 2.0 * wc * k + w0 * w0;
    next.kp = kp;
    next.b0 = 2.0 * kr * wc * k / d0;
    next.a1 = 2.0 * (w0 * w0 - k * k) / d0;
    next.a2 = (k * k - 2.0 * wc * k + w0 * w0) / d0;

    /*
     * 二阶分母 z^2+a1*z+a2 的 Jury 条件。
     * 检查极端参数溢出/舍入后的支路极点，不把它当闭环稳定性检查。
     */
    if (!isfinite(k) || !isfinite(d0) || !isfinite(next.b0) ||
        !isfinite(next.a1) || !isfinite(next.a2) ||
        fabs(next.a2) >= 1.0 ||
        1.0 + next.a1 + next.a2 <= 0.0 ||
        1.0 - next.a1 + next.a2 <= 0.0)
        return 0;
    next.ready = 1;
    *controller = next;
    return 1;
}

/* 成功返回 1；遇到非法输入/溢出返回 0，保留历史状态和输出参数。 */
static int pr_update(QuasiPr *controller, double error, double *output)
{
    double resonant, command;
    if (controller == NULL || output == NULL ||
        !controller->ready || !isfinite(error))
        return 0;

    /*
     * v[k] = b0*(e[k]-e[k-2]) - a1*v[k-1] - a2*v[k-2]。
     * 必须先用旧状态计算，再移动历史量；否则会把 e[k-2] 覆盖。
     * 分子没有 e[k-1] 项，但仍需存 e1，才能在下次形成正确 e2。
     */
    resonant = controller->b0 * (error - controller->e2) -
               controller->a1 * controller->v1 -
               controller->a2 * controller->v2;
    command = controller->kp * error + resonant;
    if (!isfinite(resonant) || !isfinite(command))
        return 0;

    controller->e2 = controller->e1;
    controller->e1 = error;
    controller->v2 = controller->v1;
    controller->v1 = resonant;
    *output = command;
    return 1;
}

/*
 * 以下复数只用于独立频率响应检查，不参与实时递推。
 * q = exp(-j*w*dt)，相当于正弦稳态下的一次采样延迟。
 */
static double complex discrete_response(const QuasiPr *controller,
                                        double w, double dt)
{
    double complex q = cexp(-I * w * dt);
    return controller->kp + controller->b0 * (1.0 - q * q) /
           (1.0 + controller->a1 * q + controller->a2 * q * q);
}

/* 按连续公式计算响应，作为离散系数验证的另一条计算路径。 */
static double complex continuous_response(double kp, double kr,
                                          double w0, double wc, double w)
{
    return kp + 2.0 * kr * wc * I * w /
           (w0 * w0 - w * w + 2.0 * wc * I * w);
}

static void check_coefficients(void)
{
    const double kp = 2.0, kr = 30.0;
    const double w0 = 2.0 * PI * 50.0, wc = 2.0 * PI * 2.0;
    const double dt = 0.0002;
    QuasiPr controller = {0};
    double output, b0, a1, a2, saved_v1;
    int i;
    assert(pr_configure(&controller, kp, kr, w0, wc, dt));
    b0 = controller.b0; a1 = controller.a1; a2 = controller.a2;

    /* 中心谐振增益是 kr；直流和 Nyquist 处谐振增益为 0。 */
    assert(cabs(discrete_response(&controller, w0, dt) - (kp + kr)) < 1e-9);
    assert(cabs(discrete_response(&controller, 0.0, dt) - kp) < 1e-12);
    assert(cabs(discrete_response(&controller, PI / dt, dt) - kp) < 1e-12);

    /*
     * 离中心频率的检查：Tustin 把数字频率 w 映射到
     * w_analog = k*tan(w*dt/2)，不能直接拿同一 w 的连续响应硬比。
     */
    for (i = 0; i < 3; ++i) {
        double w = w0 * (0.7 + 0.3 * i);
        double k = w0 / tan(w0 * dt / 2.0);
        double mapped_w = k * tan(w * dt / 2.0);
        assert(cabs(discrete_response(&controller, w, dt) -
               continuous_response(kp, kr, w0, wc, mapped_w)) < 1e-9);
    }

    /* 单位脉冲的前三拍手算，检查递推的符号和历史量移动顺序。 */
    assert(pr_update(&controller, 1.0, &output));
    assert(fabs(output - (kp + b0)) < 1e-12);
    assert(pr_update(&controller, 0.0, &output));
    assert(fabs(output - (-a1 * b0)) < 1e-12);
    assert(pr_update(&controller, 0.0, &output));
    assert(fabs(output - (a1 * a1 - a2 - 1.0) * b0) < 1e-12);

    /* 错误输入不能悄悄污染下一次计算，也不能改调用者的输出。 */
    saved_v1 = controller.v1;
    output = 123.0;
    assert(!pr_update(&controller, NAN, &output));
    assert(output == 123.0 && controller.v1 == saved_v1);
    assert(!pr_configure(&controller, kp, kr, w0, wc, 0.0));
    assert(controller.v1 == saved_v1);
    assert(!pr_configure(&controller, kp, kr, w0, wc, 0.02));
    assert(controller.v1 == saved_v1);

    /* 长时间输入常数后谐振支路衰减，说明准 PR 不提供直流积分。 */
    assert(pr_configure(&controller, kp, kr, w0, wc, dt));
    for (i = 0; i < 20000; ++i)
        assert(pr_update(&controller, 1.0, &output));
    assert(fabs(output - kp) < 1e-9);

    printf("coefficients: b0=%.12f a1=%.12f a2=%.12f\n", b0, a1, a2);
}

/*
 * 公平对比：同一 tau、采样周期、初始状态、目标幅度和时间窗口，
 * P 与准 PR 都用 kp=2，仅谐振支路的 kr 从 0 改成 30。
 * 控制输出在一个采样间隔内保持不变，故对象用精确 ZOH 更新。
 */
static Result simulate(double reference_hz, double kr)
{
    const double dt = 0.0002, tau = 0.02;
    const double kp = 2.0, w0 = 2.0 * PI * 50.0;
    const double wc = 2.0 * PI * 2.0;
    const double a = exp(-dt / tau), b = 1.0 - a;
    const double w = 2.0 * PI * reference_hz;
    QuasiPr controller = {0};
    Result result = {0};
    double y = 0.0, sum_error2 = 0.0;
    double complex z, plant, response;
    int i, count = 0;
    assert(pr_configure(&controller, kp, kr, w0, wc, dt));

    /*
     * y[k+1]=a*y[k]+b*u[k] => P(z)=b/(z-a)。
     * 闭环误差传递为 1/(1+P*C)，单位正弦的 RMS 是 1/sqrt(2)。
     * 先独立预测，再用时域结果交叉验证，避免只看“误差变小”。
     */
    z = cexp(I * w * dt);
    plant = b / (z - a);
    response = discrete_response(&controller, w, dt);
    result.predicted_rms = 1.0 / (sqrt(2.0) * cabs(1.0 + plant * response));

    for (i = 0; i < 20000; ++i) {
        double reference = sin(w * (i * dt));
        double error = reference - y;
        double output = 0.0;
        assert(pr_update(&controller, error, &output));
        result.peak_output = fmax(result.peak_output, fabs(output));

        /*
         * 前两秒只让瞬态衰减，后两秒含 50/55 Hz 的整数个周期。
         * error 属于时刻 k，所以必须先统计，再更新到 y[k+1]。
         */
        if (i >= 10000) {
            sum_error2 += error * error;
            ++count;
        }
        y = a * y + b * output;
        assert(isfinite(y));
    }
    result.rms_error = sqrt(sum_error2 / count);
    assert(fabs(result.rms_error - result.predicted_rms) < 1e-8);
    return result;
}

int main(void)
{
    Result p50, pr50, p55, pr55;
    check_coefficients();
    p50 = simulate(50.0, 0.0);
    pr50 = simulate(50.0, 30.0);
    p55 = simulate(55.0, 0.0);
    pr55 = simulate(55.0, 30.0);

    /* 验证本组参数的改善及失谐退化，不把它推广为任意对象的保证。 */
    assert(pr50.rms_error < 0.25 * p50.rms_error);
    assert(pr55.rms_error > 4.0 * pr50.rms_error);
    assert(pr55.rms_error < p55.rms_error);
    assert(pr50.peak_output < 8.0 && pr55.peak_output < 12.0);

    printf("50 Hz P:  error_rms=%.9f predicted=%.9f peak_u=%.6f\n",
           p50.rms_error, p50.predicted_rms, p50.peak_output);
    printf("50 Hz PR: error_rms=%.9f predicted=%.9f peak_u=%.6f\n",
           pr50.rms_error, pr50.predicted_rms, pr50.peak_output);
    printf("55 Hz P:  error_rms=%.9f predicted=%.9f peak_u=%.6f\n",
           p55.rms_error, p55.predicted_rms, p55.peak_output);
    printf("55 Hz PR: error_rms=%.9f predicted=%.9f peak_u=%.6f\n",
           pr55.rms_error, pr55.predicted_rms, pr55.peak_output);
    puts("PASS: quasi-PR coefficients, recurrence, frequency response and tracking.");
    return 0;
}

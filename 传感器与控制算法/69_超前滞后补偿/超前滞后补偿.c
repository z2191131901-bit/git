/*
 * 69 超前滞后补偿：把连续传递函数变成每个周期调用一次的 C 函数。
 * 标准 C99。所有控制信号按归一化值处理，时间单位是秒。
 * 模型仿真仅说明本组参数下的行为，没有连接实际执行器。
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>

#define PI 3.14159265358979323846

/* C(z)=(b0+b1*z^-1)/(1+a1*z^-1)，历史量按一次采样保存。 */
typedef struct {
    double b0, b1, a1;
    double previous_input, previous_output;
} Compensator;

typedef struct {
    double real, imaginary;
} Response;

typedef struct {
    double final_output;
    double time_to_90;       /* 首次达到自身最终值 90% 的时间，单位 s。 */
    double settling_time;   /* 最后离开自身最终值 ±2% 带之后的采样时刻。 */
    double maximum_input;   /* 控制输入的绝对值峰值。 */
} Simulation;

/*
 * 初始化 C(s)=gain*(1+zero_time*s)/(1+pole_time*s)。
 * zero_time、pole_time、dt 单位都是 s；gain 是输出/输入的比例。
 * zero_time>pole_time：超前；zero_time<pole_time：滞后。
 * 使用等价的 dt+2*T 形式，避免直接计算 2/dt。
 * 失败返回 0，并保持原对象；成功清零历史量。
 */
static int compensator_init(Compensator *filter, double gain,
                            double zero_time, double pole_time, double dt)
{
    Compensator next;
    double denominator;
    if (filter == NULL || !isfinite(gain) || gain < 0.0 ||
        !isfinite(zero_time) || zero_time < 0.0 ||
        !isfinite(pole_time) || pole_time <= 0.0 ||
        !isfinite(dt) || dt <= 0.0) {
        return 0;
    }
    denominator = dt + 2.0 * pole_time;
    next.b0 = gain * (dt + 2.0 * zero_time) / denominator;
    next.b1 = gain * (dt - 2.0 * zero_time) / denominator;
    next.a1 = (dt - 2.0 * pole_time) / denominator;
    next.previous_input = 0.0;
    next.previous_output = 0.0;
    /*
     * 理论正时间常数对应 |a1|<1。
     * 极端量级可能因浮点舍入得到单位圆上的极点，此处也拒绝。
     */
    if (!isfinite(denominator) || !isfinite(next.b0) ||
        !isfinite(next.b1) || !isfinite(next.a1) ||
        fabs(next.a1) >= 1.0) {
        return 0;
    }
    *filter = next;
    return 1;
}

/*
 * 先算新输出，再更新历史值。千万不要先覆盖 previous_input。
 * output 不能指向 filter 内部成员；调用者需先成功初始化 filter。
 * 非有限输入或运算溢出均返回 0，不改历史量，也不改 output。
 */
static int compensator_step(Compensator *filter, double input, double *output)
{
    double next_output;
    if (filter == NULL || output == NULL || !isfinite(input)) {
        return 0;
    }
    next_output = filter->b0 * input +
                  filter->b1 * filter->previous_input -
                  filter->a1 * filter->previous_output;
    if (!isfinite(next_output)) {
        return 0;
    }
    filter->previous_input = input;
    filter->previous_output = next_output;
    *output = next_output;
    return 1;
}

/* 检查用：在 z=exp(j*theta) 处计算频率响应，theta 单位 rad/sample。 */
static Response discrete_response(const Compensator *filter, double theta)
{
    const double nr = filter->b0 + filter->b1 * cos(theta);
    const double ni = -filter->b1 * sin(theta);
    const double dr = 1.0 + filter->a1 * cos(theta);
    const double di = -filter->a1 * sin(theta);
    const double divisor = dr * dr + di * di;
    Response answer = {(nr * dr + ni * di) / divisor,
                       (ni * dr - nr * di) / divisor};
    return answer;
}

/* 检查用：连续公式 K*(1+j*w*Tz)/(1+j*w*Tp)。 */
static Response continuous_response(double gain, double tz, double tp, double w)
{
    const double divisor = 1.0 + w * w * tp * tp;
    Response answer = {gain * (1.0 + w * w * tz * tp) / divisor,
                       gain * w * (tz - tp) / divisor};
    return answer;
}

static void check_frequency_response(void)
{
    Compensator lead, lag;
    const double dt = 0.01;
    const double omega = 1.0 / sqrt(0.5 * 0.1);
    const double theta = 2.0 * atan(omega * dt / 2.0);
    Response digital, analog, lag_middle, dc, nyquist;
    double cosine_coefficient = 0.0, sine_coefficient = 0.0;
    const int period = 200;
    const int samples = 20 * period;
    int index;

    assert(compensator_init(&lead, 2.0, 0.5, 0.1, dt));
    assert(compensator_init(&lag, 5.0, 0.5, 2.5, dt));
    dc = discrete_response(&lead, 0.0);
    nyquist = discrete_response(&lead, PI);
    assert(fabs(dc.real - 2.0) < 1e-12);
    assert(fabs(nyquist.real - 10.0) < 1e-12);
    assert(fabs(discrete_response(&lag, 0.0).real - 5.0) < 1e-11);
    assert(fabs(discrete_response(&lag, PI).real - 1.0) < 1e-12);

    /* 按 Tustin 频率映射比较，而非误把模拟 w*dt 直接当成 theta。 */
    digital = discrete_response(&lead, theta);
    analog = continuous_response(2.0, 0.5, 0.1, omega);
    assert(fabs(digital.real - analog.real) < 1e-12);
    assert(fabs(digital.imaginary - analog.imaginary) < 1e-12);
    assert(atan2(digital.imaginary, digital.real) > 0.0);
    lag_middle = discrete_response(&lag,
        2.0 * atan(dt / (2.0 * sqrt(0.5 * 2.5))));
    assert(atan2(lag_middle.imaginary, lag_middle.real) < 0.0);

    /*
     * 用实际递推处理正弦输入，再做整周期投影。
     * 对输入 cos(theta*k)，输出为 Re(H)*cos - Im(H)*sin。
     * 先运行 10 个周期消除启动暂态，之后检查振幅和相位两部分。
     */
    for (index = -10 * period; index < samples; ++index) {
        const double phase = 2.0 * PI * index / period;
        double output;
        assert(compensator_step(&lead, cos(phase), &output));
        if (index >= 0) {
            cosine_coefficient += output * cos(phase);
            sine_coefficient += output * sin(phase);
        }
    }
    cosine_coefficient *= 2.0 / samples;
    sine_coefficient *= 2.0 / samples;
    digital = discrete_response(&lead, 2.0 * PI / period);
    assert(fabs(cosine_coefficient - digital.real) < 1e-10);
    assert(fabs(sine_coefficient + digital.imaginary) < 1e-10);
    printf("frequency: lead maximum phase = %.6f deg, lag middle phase = %.6f deg\n",
           atan2(analog.imaginary, analog.real) * 180.0 / PI,
           atan2(lag_middle.imaginary, lag_middle.real) * 180.0 / PI);
}

static void check_boundaries(void)
{
    Compensator filter, snapshot;
    double output = 123.0;
    int index;
    assert(compensator_init(&filter, 2.0, 0.5, 0.1, 0.01));
    /* 零历史下，单位阶跃首输出就是 b0=2*1.01/0.21。 */
    assert(compensator_step(&filter, 1.0, &output));
    assert(fabs(output - 9.619047619047619) < 1e-12);
    snapshot = filter;
    assert(!compensator_step(&filter, NAN, &output));
    assert(output == snapshot.previous_output);
    assert(filter.previous_input == snapshot.previous_input);
    assert(filter.previous_output == snapshot.previous_output);
    assert(!compensator_init(&filter, 1.0, 0.1, 0.0, 0.01));
    assert(!compensator_init(&filter, 1.0, 0.1, 0.2, 0.0));
    assert(!compensator_init(&filter, -1.0, 0.1, 0.2, 0.01));
    assert(filter.b0 == snapshot.b0 && filter.a1 == snapshot.a1);
    assert(filter.previous_output == snapshot.previous_output);

    /* Tz=Tp 时，零极点抵消，应对任意输入退化为纯比例 gain。 */
    assert(compensator_init(&filter, 3.0, 0.2, 0.2, 0.01));
    for (index = 0; index < 200; ++index) {
        const double input = sin(index * 0.13) + (index % 3);
        assert(compensator_step(&filter, input, &output));
        assert(fabs(output - 3.0 * input) < 1e-12);
    }
    /* 零增益合法；零点时间为零时也可退化为一阶低通。 */
    assert(compensator_init(&filter, 0.0, 0.5, 0.1, 0.01));
    assert(compensator_step(&filter, 1.0, &output) && output == 0.0);
    assert(compensator_init(&filter, 1.0, 0.0, 0.1, 0.01));
    for (index = 0; index < 1000; ++index) {
        assert(compensator_step(&filter, 1.0, &output));
    }
    assert(fabs(output - 1.0) < 1e-12);
}

/*
 * 四组使用同一对象 y'=-y+u、同一个单位阶跃、相同 dt 与初值。
 * 对象每个周期内输入保持不变，因此用解析离散式推进。
 * use_lead=0：第一段是比例 2；=1：第一段是直流增益 2 的超前。
 * use_lag=1：后接直流增益 5、高频增益 1 的滞后。
 * 本演示不施加限幅，便于检查线性理论；会明确报告控制量峰值。
 */
static Simulation simulate(int use_lead, int use_lag)
{
    const double dt = 0.01;
    const double plant_a = exp(-dt);
    const double dc_gain = use_lag ? 10.0 : 2.0;
    const double expected_output = dc_gain / (1.0 + dc_gain);
    Compensator lead, lag;
    Simulation result = {0.0, -1.0, 0.0, 0.0};
    double y = 0.0;
    int index;
    assert(compensator_init(&lead, 2.0, 0.5, 0.1, dt));
    assert(compensator_init(&lag, 5.0, 0.5, 2.5, dt));
    for (index = 0; index < 6000; ++index) {
        const double error = 1.0 - y;
        double first = 2.0 * error;
        double input;
        if (use_lead) {
            assert(compensator_step(&lead, error, &first));
        }
        input = first;
        if (use_lag) {
            assert(compensator_step(&lag, first, &input));
        }
        if (fabs(input) > result.maximum_input) {
            result.maximum_input = fabs(input);
        }
        y = plant_a * y + (1.0 - plant_a) * input;
        assert(isfinite(y) && fabs(y) < 5.0);
        if (result.time_to_90 < 0.0 && y >= 0.9 * expected_output) {
            result.time_to_90 = (index + 1) * dt;
        }
        if (fabs(y - expected_output) > 0.02 * fabs(expected_output)) {
            /* 下一个采样时刻可能开始进入并永久保持在误差带内。 */
            result.settling_time = (index + 2) * dt;
        }
    }
    result.final_output = y;
    assert(fabs(y - expected_output) < 1e-8);
    assert(result.time_to_90 > 0.0 && result.settling_time < 60.0);
    return result;
}

int main(void)
{
    const char *names[4] = {"P", "lead", "P+lag", "lead+lag"};
    Simulation result[4];
    int index;
    check_frequency_response();
    check_boundaries();
    result[0] = simulate(0, 0);
    result[1] = simulate(1, 0);
    result[2] = simulate(0, 1);
    result[3] = simulate(1, 1);
    for (index = 0; index < 4; ++index) {
        printf("%-8s final=%.9f error=%.9f t90=%.3f s settle=%.3f s max|u|=%.6f\n",
               names[index], result[index].final_output,
               1.0 - result[index].final_output, result[index].time_to_90,
               result[index].settling_time, result[index].maximum_input);
    }
    /* 只断言当前模型与这组参数的结论，不声称对任意对象均成立。 */
    assert(fabs(result[0].final_output - result[1].final_output) < 1e-8);
    assert(result[1].time_to_90 < result[0].time_to_90);
    assert(result[1].maximum_input > result[0].maximum_input);
    assert(1.0 - result[2].final_output < 1.0 - result[0].final_output);
    assert(1.0 - result[3].final_output < 1.0 - result[1].final_output);
    puts("PASS: lead/lag recurrence, frequency response, boundaries and simulations.");
    return 0;
}

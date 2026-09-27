/*
 * 第72篇：ZV / ZVD 输入整形。
 * 整形的是给对象的目标指令，不是对测量做滤波。
 * 用两个或三个延迟副本的加权和，减小已知振动模态的残振。
 * 本例使用整数拍延迟，讲解同时说明延迟量化及频率失配的影响。
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>

#define PI 3.14159265358979323846
#define HISTORY 2048

typedef struct {
    double weight[3];       /* 脉冲权重，和为1，保证恒定输入最终不变。 */
    int delay[3];           /* 各脉冲延迟，单位为采样拍数。 */
    int count;              /* 1=原指令，2=ZV，3=ZVD。 */
    int write;              /* 下一次写入当前输入的环形数组位置。 */
    double history[HISTORY]; /* 启动前的输入假定为0。 */
} Shaper;

/*
 * kind=0:不整形；kind=1:ZV；kind=2:ZVD。
 * fn是设计模态的无阻尼自然频率Hz，zeta为阻尼比，dt为采样秒数。
 * 仅支持欠阻尼模态0<=zeta<1。初始化成功会清空历史。
 */
static int shaper_init(Shaper *s, int kind, double fn, double zeta, double dt)
{
    Shaper next = {0};
    double half_period_samples, decay, denominator;
    int n;
    if (!s || kind < 0 || kind > 2 || !isfinite(fn) || fn <= 0.0 ||
        !isfinite(zeta) || zeta < 0.0 || zeta >= 1.0 ||
        !isfinite(dt) || dt <= 0.0)
        return 0;
    half_period_samples = 1.0 / (2.0 * fn * sqrt(1.0 - zeta * zeta) * dt);
    /*
     * 四舍五入到整数拍。保留完整ZVD所需的2*n拍，即便当前选ZV。
     * 先限制范围再转int，避免过大浮点数转换为整数。
     */
    if (!isfinite(half_period_samples) || half_period_samples < 0.5 ||
        half_period_samples >= (HISTORY - 1) / 2.0)
        return 0;
    n = (int)floor(half_period_samples + 0.5);
    decay = exp(-zeta * PI / sqrt(1.0 - zeta * zeta));
    denominator = 1.0 + decay;
    next.count = kind + 1;
    next.weight[0] = 1.0;
    if (kind == 1) {
        next.weight[0] = 1.0 / denominator;
        next.weight[1] = decay / denominator;
        next.delay[1] = n;
    } else if (kind == 2) {
        /* ZVD相当于两个相同ZV整形器卷积，系数遵循平方展开。 */
        next.weight[0] = 1.0 / (denominator * denominator);
        next.weight[1] = 2.0 * decay / (denominator * denominator);
        next.weight[2] = decay * decay / (denominator * denominator);
        next.delay[1] = n;
        next.delay[2] = 2 * n;
    }
    *s = next;
    return 1;
}

/* s必须由init成功配置，外部不应改写权重、延迟和数组索引。 */
static int shaper_step(Shaper *s, double input, double *output)
{
    double sum = 0.0;
    if (!s || !output || s->count < 1 || s->count > 3 || !isfinite(input))
        return 0;
    for (int j = 0; j < s->count; ++j) {
        int index = (s->write + HISTORY - s->delay[j]) % HISTORY;
        /* 零延迟项使用当前输入，其余项读取尚未覆盖的旧输入。 */
        double sample = s->delay[j] == 0 ? input : s->history[index];
        sum += s->weight[j] * sample;
    }
    if (!isfinite(sum))
        return 0;
    s->history[s->write] = input;
    s->write = (s->write + 1) % HISTORY;
    *output = sum;
    return 1;
}

/*
 * 二阶对象：y''+2*zeta*wn*y'+wn^2*y=wn^2*u。
 * 每拍u恒定，下面是精确零阶保持更新，不依赖数值积分步长近似。
 */
static void plant_step(double *y, double *v, double u,
                       double wn, double zeta, double dt)
{
    double wd = wn * sqrt(1.0 - zeta * zeta);
    double decay = exp(-zeta * wn * dt);
    double cs = cos(wd * dt), sn = sin(wd * dt);
    double error = *y - u;
    double next_y = u + decay * ((cs + zeta * wn / wd * sn) * error + sn / wd * *v);
    double next_v = decay * (-wn * wn / wd * sn * error +
                            (cs - zeta * wn / wd * sn) * *v);
    *y = next_y;
    *v = next_v;
}

/* 连续单位阶跃解析解，用另一条路径核对采样仿真的y和v。 */
static void step_response(double t, double wn, double zeta, double *y, double *v)
{
    double wd = wn * sqrt(1.0 - zeta * zeta);
    double decay = exp(-zeta * wn * t);
    *y = 1.0 - decay * (cos(wd * t) + zeta * wn / wd * sin(wd * t));
    *v = decay * wn * wn / wd * sin(wd * t);
}

static double residual(int kind, double actual_hz)
{
    const double dt = 0.001, zeta = 0.05, wn = 2.0 * PI * actual_hz;
    Shaper s, full;
    double y = 0.0, v = 0.0, exact_y = 0.0, exact_v = 0.0;
    assert(shaper_init(&s, kind, 5.0, zeta, dt));
    assert(shaper_init(&full, 2, 5.0, zeta, dt));
    /* 所有配置都在同一时刻比较：最迟ZVD脉冲生效后再走一拍。 */
    int steps = full.delay[2] + 1;
    for (int k = 0; k < steps; ++k) {
        double command;
        assert(shaper_step(&s, 1.0, &command));
        plant_step(&y, &v, command, wn, zeta, dt);
    }
    for (int j = 0; j < s.count; ++j) {
        double component_y, component_v;
        double elapsed = (steps - s.delay[j]) * dt;
        step_response(elapsed, wn, zeta, &component_y, &component_v);
        exact_y += s.weight[j] * component_y;
        exact_v += s.weight[j] * component_v;
    }
    assert(fabs(y - exact_y) < 1e-11 && fabs(v - exact_v) < 1e-10);

    /*
     * 最终指令已经等于1，后续误差是衰减正弦。
     * 用位置和速度一起恢复当前包络幅值，避免刚好取到过零点。
     */
    double wd = wn * sqrt(1.0 - zeta * zeta);
    return hypot(y - 1.0, (v + zeta * wn * (y - 1.0)) / wd);
}

int main(void)
{
    Shaper s;
    const double zeta = 0.05;
    const double decay = exp(-zeta * PI / sqrt(1.0 - zeta * zeta));
    double output, sum = 0.0;
    assert(shaper_init(&s, 2, 5.0, zeta, 0.001));
    assert(s.delay[1] == 100 && s.delay[2] == 200);
    assert(fabs(s.weight[0] - 1.0 / ((1.0 + decay) * (1.0 + decay))) < 1e-14);

    /* 单位脉冲应当只在0、100、200拍出现，且三个权重加起来为1。 */
    for (int k = 0; k <= 200; ++k) {
        assert(shaper_step(&s, k == 0 ? 1.0 : 0.0, &output));
        double expected = k == 0 ? s.weight[0] :
                          k == 100 ? s.weight[1] : k == 200 ? s.weight[2] : 0.0;
        assert(fabs(output - expected) < 1e-14);
        sum += output;
    }
    assert(fabs(sum - 1.0) < 1e-14);

    /* 绕过数组末端多次，验证延迟索引和恒定指令的直流增益。 */
    for (int k = 0; k < 3 * HISTORY; ++k)
        assert(shaper_step(&s, 1.0, &output));
    assert(fabs(output - 1.0) < 1e-14);
    int saved_write = s.write;
    output = 123.0;
    assert(!shaper_step(&s, NAN, &output));
    assert(s.write == saved_write && output == 123.0);
    assert(!shaper_init(&s, 2, 0.0, zeta, 0.001));
    assert(s.write == saved_write);
    assert(!shaper_init(&s, 2, 5.0, 1.0, 0.001));
    assert(!shaper_init(&s, 2, 0.01, zeta, 0.001));

    double raw = residual(0, 5.0);
    double zv = residual(1, 5.0);
    double zvd = residual(2, 5.0);
    double raw_offset = residual(0, 5.5);
    double zv_offset = residual(1, 5.5);
    double zvd_offset = residual(2, 5.5);
    assert(zv < 0.005 * raw && zvd < 0.00003 * raw);
    assert(zvd_offset < 0.25 * zv_offset && zv_offset < raw_offset);
    printf("shaper: common_time=0.201 s; delay ZV=0.100 s ZVD=0.200 s\n");
    printf("shaper: 5 Hz residual raw=%.9f ZV=%.9f ZVD=%.9f\n", raw, zv, zvd);
    printf("shaper: 5.5 Hz residual raw=%.9f ZV=%.9f ZVD=%.9f\n",
           raw_offset, zv_offset, zvd_offset);
    puts("PASS: impulse weights, ring wrap, analytical response and frequency mismatch.");
    return 0;
}

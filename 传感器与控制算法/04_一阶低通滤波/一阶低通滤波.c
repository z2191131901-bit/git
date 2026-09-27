/*
 * 一阶低通：输出 += alpha * (输入 - 输出)。
 * tau 是时间常数（秒），dt 是采样间隔（秒）。
 * alpha = dt / (tau + dt)，采用后向欧拉离散化。
 */
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

typedef struct
{
    float value;     /* 上一次滤波结果；下一次更新仍会用到。 */
    int initialized; /* 0：尚无有效样本，1：已初始化。 */
} LowPass;

/* 使用 {0} 初始化；输出保存在 state->value。失败时不修改状态。 */
static int low_pass_update(LowPass *state, float input, float tau, float dt)
{
    if (state == NULL || !isfinite(input) || !isfinite(tau) ||
        !isfinite(dt) || tau <= 0 || dt <= 0)
        return 0;

    /*
     * 第一个有效样本直接作为起点。
     * 否则温度本来是 25 度，却会从默认的 0 度缓慢爬升。
     */
    if (!state->initialized)
    {
        state->value = input;
        state->initialized = 1;
        return 1;
    }

    /* double 中间量降低 tau+dt 和 input-value 的 float 溢出风险。 */
    const double alpha = (double)dt / ((double)tau + dt);
    const double next = state->value + alpha * ((double)input - state->value);
    if (!isfinite(next))
        return 0;
    state->value = (float)next;
    return 1;
}

int main(void)
{
    LowPass filter = {0};
    int ok = low_pass_update(&filter, 10, 0.03f, 0.01f);
    assert(ok && filter.value == 10); /* 首样本直接输出。 */

    ok = low_pass_update(&filter, 14, 0.03f, 0.01f);
    assert(ok && fabsf(filter.value - 11) < 1e-5f);
    printf("first step=%.6f (expected 11)\n", filter.value);

    ok = low_pass_update(&filter, 14, 0.03f, 0.01f);
    assert(ok && fabsf(filter.value - 11.75f) < 1e-5f);

    const float saved = filter.value;
    ok = low_pass_update(&filter, 14, 0.03f, 0);
    assert(!ok && filter.value == saved);

    /* 恒定输入下应逐渐收敛到输入值。 */
    for (int i = 0; i < 100; ++i)
    {
        ok = low_pass_update(&filter, 14, 0.03f, 0.01f);
        assert(ok);
    }
    assert(fabsf(filter.value - 14) < 1e-4f);
    puts("low_pass: PASS");
    return 0;
}

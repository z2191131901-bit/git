/*
 * 滑动平均：保存最近 4 个有效采样，输出它们的平均值。
 * 读法：先看 main 的输入，再看结构体，最后看 update 的四个步骤。
 * 输入、输出单位相同；每个传感器通道需要自己的一份状态。
 */
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

#define WINDOW_SIZE 4

typedef struct
{
    float samples[WINDOW_SIZE]; /* 历史读数。下标范围是 0 到 3。 */
    size_t next;                /* 下一次写入的位置；窗口满后也是最旧值的位置。 */
    size_t count;               /* 已经收到多少个有效样本，最大为 4。 */
    double sum;                 /* 用 double 累加，减小长期加减的舍入误差。 */
} MovingAverage;

/*
 * state 指向跨采样保存的状态，output 指向存放结果的 float。
 * 首次使用必须用 {0} 初始化 state；output 不得指向 state 内部。
 * 返回 1 为成功，0 为失败；失败时状态和 output 都保持不变。
 */
static int moving_average_update(MovingAverage *state, float input, float *output)
{
    if (state == NULL || output == NULL || !isfinite(input))
        return 0;

    /* ① 窗口已满时，旧值即将被覆盖，必须先从总和里扣掉。 */
    if (state->count == WINDOW_SIZE)
        state->sum -= state->samples[state->next];
    else
        state->count++; /* 尚未填满时，只增加有效数据数量。 */

    /* ② 新值加入窗口和总和。减旧加新避免每次遍历整个数组。 */
    state->samples[state->next] = input;
    state->sum += input;

    /*
     * ③ 下标加 1，并在到达 4 时回到 0。
     * % 是取余运算，因此 next 依次为 0、1、2、3、0……
     */
    state->next = (state->next + 1) % WINDOW_SIZE;

    /* ④ 启动阶段除以有效数量，而不是除以固定的窗口长度。 */
    *output = (float)(state->sum / (double)state->count);
    return 1;
}

int main(void)
{
    MovingAverage filter = {0}; /* 只初始化一次，不能放在采样循环里。 */
    const float input[] = {10, 12, 14, 16, 18, 20};
    const float expected[] = {10, 11, 12, 13, 15, 17};
    float output = 0;

    for (size_t i = 0; i < sizeof(input) / sizeof(input[0]); ++i)
    {
        int ok = moving_average_update(&filter, input[i], &output);
        assert(ok);
        assert(fabsf(output - expected[i]) < 1e-6f);
        printf("input=%.1f average=%.1f\n", input[i], output);
    }

    /* 非法读数不应占用窗口；后续有效读数仍接在原序列后。 */
    int ok = moving_average_update(&filter, NAN, &output);
    assert(!ok && output == 17);
    ok = moving_average_update(&filter, 22, &output);
    assert(ok && output == 19);
    puts("moving_average: PASS");
    return 0;
}

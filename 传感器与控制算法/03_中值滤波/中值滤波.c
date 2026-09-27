/*
 * 中值滤波：对最近 5 个有效读数排序，取中间值。
 * 为了便于学习，用插入排序；不依赖第三方库，不申请堆内存。
 */
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

#define WINDOW_SIZE 5
typedef struct
{
    float samples[WINDOW_SIZE]; /* 保留采样时间顺序的循环缓冲区。 */
    size_t next;
    size_t count;
} MedianFilter;

/* state 需先用 {0} 初始化；output 不得指向 state 内部。 */
static int median_update(MedianFilter *state, float input, float *output)
{
    if (state == NULL || output == NULL || !isfinite(input))
        return 0;

    state->samples[state->next] = input;
    state->next = (state->next + 1) % WINDOW_SIZE;
    if (state->count < WINDOW_SIZE)
        state->count++;

    /*
     * 只排序副本，不能打乱历史缓冲区。
     * 原缓冲区的 next 还需要准确找到下次应替换的最旧样本。
     */
    float sorted[WINDOW_SIZE];
    for (size_t i = 0; i < state->count; ++i)
        sorted[i] = state->samples[i];

    /*
     * 插入排序：左侧已有序，从右侧拿一个数 key 放到正确位置。
     * j > 0 必须先判断，避免 size_t 无符号数在 0-1 时下溢。
     */
    for (size_t i = 1; i < state->count; ++i)
    {
        const float key = sorted[i];
        size_t j = i;
        while (j > 0 && sorted[j - 1] > key)
        {
            sorted[j] = sorted[j - 1]; /* 大数向右移，给 key 腾位置。 */
            --j;
        }
        sorted[j] = key;
    }

    const size_t middle = state->count / 2; /* 整数除法舍去小数。 */
    if (state->count % 2 == 1)
        *output = sorted[middle]; /* 奇数个：取正中间。 */
    else
        /* 启动时可能是偶数个：取中间两数平均。先转 double 避免求和溢出。 */
        *output = (float)(((double)sorted[middle - 1] + sorted[middle]) / 2.0);
    return 1;
}

int main(void)
{
    MedianFilter filter = {0};
    const float input[] = {10, 11, 100, 12, 13, 14, 15};
    const float expected[] = {10, 10.5f, 11, 11.5f, 12, 13, 14};
    float output = 0;
    for (size_t i = 0; i < sizeof(input) / sizeof(input[0]); ++i)
    {
        int ok = median_update(&filter, input[i], &output);
        assert(ok && fabsf(output - expected[i]) < 1e-6f);
        printf("input=%.1f median=%.1f\n", input[i], output);
    }
    int ok = median_update(&filter, INFINITY, &output);
    assert(!ok && output == 14);
    ok = median_update(&filter, 16, &output);
    assert(ok && output == 14); /* 窗口此时为 12、13、14、15、16。 */
    puts("median: PASS");
    return 0;
}

/*
 * 静止零偏标定：已确认静止时，分别求三轴读数的平均值。
 * 本例不自动判断静止！运动时求出的均值不能当零偏。
 * 输入、零偏及修正结果均使用 rad/s。
 */
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

/*
 * samples[][3]：每行一个采样，每行固定包含 x、y、z 三个分量。
 * count：行数。bias[3]：输出三轴零偏，不能和输入数组重叠。
 * 返回 1 成功；0 表示样本不足或数据无效，输出保持不变。
 */
static int gyro_bias_estimate(const float samples[][3], size_t count, float bias[3])
{
    if (samples == NULL || bias == NULL || count < 2)
        return 0;

    /* 先求所有样本的总和，使用 double 保存累加结果。 */
    double sum[3] = {0, 0, 0};
    for (size_t i = 0; i < count; ++i)
    {
        for (size_t axis = 0; axis < 3; ++axis)
        {
            if (!isfinite(samples[i][axis]))
                return 0;
            sum[axis] += samples[i][axis];
        }
    }

    /* 先准备完整结果，再写输出，避免某一轴失败时留下半份结果。 */
    float next[3];
    for (size_t axis = 0; axis < 3; ++axis)
    {
        const double mean = sum[axis] / (double)count;
        if (!isfinite(mean))
            return 0;
        next[axis] = (float)mean;
    }
    for (size_t axis = 0; axis < 3; ++axis)
        bias[axis] = next[axis];
    return 1;
}

int main(void)
{
    /* 四个静止样本只是为了手算，真实标定一般采集更多样本。 */
    const float samples[][3] = {
        {0.01f, -0.02f, 0.03f},
        {0.02f, -0.01f, 0.04f},
        {0.03f, -0.02f, 0.05f},
        {0.02f, -0.03f, 0.04f}
    };
    float bias[3] = {0, 0, 0};
    int ok = gyro_bias_estimate(samples, 4, bias);
    assert(ok);
    assert(fabsf(bias[0] - 0.02f) < 1e-6f);
    assert(fabsf(bias[1] + 0.02f) < 1e-6f);
    assert(fabsf(bias[2] - 0.04f) < 1e-6f);

    /*
     * 新读数 = 真实角速度 + 零偏 + 噪声。
     * 去零偏就是减去估计偏差，不能把噪声也完全消掉。
     */
    const float measured_z = 0.14f;
    const float corrected_z = measured_z - bias[2];
    assert(fabsf(corrected_z - 0.10f) < 1e-6f);
    printf("bias=(%.4f, %.4f, %.4f), corrected_z=%.4f rad/s\n",
           bias[0], bias[1], bias[2], corrected_z);

    const float bad[][3] = {{0, 0, 0}, {0, NAN, 0}};
    ok = gyro_bias_estimate(bad, 2, bias);
    assert(!ok && fabsf(bias[2] - 0.04f) < 1e-6f);
    ok = gyro_bias_estimate(samples, 1, bias);
    assert(!ok);
    puts("gyro_bias: PASS");
    return 0;
}

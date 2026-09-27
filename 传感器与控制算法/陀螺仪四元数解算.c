/*
 * 入门样稿：用三轴陀螺仪更新四元数。
 * 配套讲解：01_陀螺仪四元数解算.md
 *
 * 约定：右手坐标系；四元数顺序为 w、x、y、z；使用 Hamilton 乘法。
 * 四元数表示从机体坐标系到世界坐标系的旋转。
 * 输入角速度在机体坐标系中表达，单位 rad/s；时间单位 s。
 * 本例只有陀螺仪积分，不包含加速度计校正，也不能消除漂移。
 */
#include <math.h>
#include <stddef.h>
#include <stdio.h>

typedef struct
{
    float w;  /* 标量部分；初始姿态取单位四元数时为 1。 */
    float x;  /* 以下三个数是向量部分，不是三个欧拉角。 */
    float y;
    float z;
} Quaternion;

/* 将常见的“度/秒”换算为算法需要的“弧度/秒”。 */
static float degrees_to_radians(float degrees)
{
    const float pi = 3.14159265358979323846f;
    return degrees * pi / 180.0f;
}

/*
 * 更新一次姿态。
 * q：输入上一次姿态，成功后保存本次姿态；首次使用可设为 {1,0,0,0}。
 * gyro_x/y/z：已经换算好单位、对齐坐标轴的三轴角速度，单位 rad/s。
 * dt：两次采样之间的实际时间间隔，单位 s。
 * 返回 1 表示成功；返回 0 表示输入或计算无效，此时不修改 q。
 */
static int quaternion_update(Quaternion *q,
                             float gyro_x, float gyro_y, float gyro_z,
                             float dt)
{
    if (q == NULL || !isfinite(dt) || dt <= 0.0f ||
        !isfinite(gyro_x) || !isfinite(gyro_y) || !isfinite(gyro_z))
    {
        return 0;
    }

    /* 第一步：保存旧值。四个新分量必须由同一时刻的旧姿态计算。 */
    const float old_w = q->w;
    const float old_x = q->x;
    const float old_y = q->y;
    const float old_z = q->z;

    /* 第二步：计算四元数变化率，对应 q_dot = 0.5 * q ⊗ [0, ω]。 */
    const float derivative_w = 0.5f * (-old_x * gyro_x - old_y * gyro_y - old_z * gyro_z);
    const float derivative_x = 0.5f * ( old_w * gyro_x + old_y * gyro_z - old_z * gyro_y);
    const float derivative_y = 0.5f * ( old_w * gyro_y - old_x * gyro_z + old_z * gyro_x);
    const float derivative_z = 0.5f * ( old_w * gyro_z + old_x * gyro_y - old_y * gyro_x);

    /* 第三步：一阶欧拉积分。新值 ≈ 旧值 + 变化率 × 时间间隔。 */
    Quaternion next;
    next.w = old_w + derivative_w * dt;
    next.x = old_x + derivative_x * dt;
    next.y = old_y + derivative_y * dt;
    next.z = old_z + derivative_z * dt;

    /* 第四步：归一化，修正数值积分造成的长度偏差。 */
    const float norm_squared = next.w * next.w + next.x * next.x
                             + next.y * next.y + next.z * next.z;
    if (!isfinite(norm_squared) || norm_squared <= 0.0f)
    {
        return 0;
    }

    const float inverse_norm = 1.0f / sqrtf(norm_squared);
    next.w *= inverse_norm;
    next.x *= inverse_norm;
    next.y *= inverse_norm;
    next.z *= inverse_norm;

    /* 全部计算成功后，再一次性写回姿态。 */
    *q = next;
    return 1;
}

int main(void)
{
    /* 把启动时的方向定义为参考方向，不代表自动找到了地理方向。 */
    Quaternion attitude = {1.0f, 0.0f, 0.0f, 0.0f};

    /* 模拟绕 Z 轴正向匀速旋转：90 度/秒，持续 1 秒。 */
    const float gyro_z = degrees_to_radians(90.0f);
    const float sample_period = 0.01f;  /* 100 Hz 采样对应 0.01 秒。 */

    for (int sample = 0; sample < 100; ++sample)
    {
        if (!quaternion_update(&attitude, 0.0f, 0.0f, gyro_z, sample_period))
        {
            fprintf(stderr, "Quaternion update failed.\n");
            return 1;
        }
    }

    printf("w=%.6f, x=%.6f, y=%.6f, z=%.6f\n",
           attitude.w, attitude.x, attitude.y, attitude.z);

    /* 理论结果为 [cos(45°), 0, 0, sin(45°)]；允许少量积分误差。 */
    const float expected = sqrtf(0.5f);
    const float tolerance = 0.0001f;
    if (fabsf(attitude.w - expected) > tolerance ||
        fabsf(attitude.x) > tolerance || fabsf(attitude.y) > tolerance ||
        fabsf(attitude.z - expected) > tolerance)
    {
        fprintf(stderr, "Reference case failed.\n");
        return 1;
    }

    return 0;
}

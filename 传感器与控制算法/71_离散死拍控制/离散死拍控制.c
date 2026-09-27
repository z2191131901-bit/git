/*
 * 第71篇：一阶离散死拍控制。
 *
 * 已知对象：y[k+1] = a*y[k] + b*u[k]。
 * 希望下一拍到达 r[k+1]，解出 u[k] = (r[k+1] - a*y[k]) / b。
 * 模型准确、没有额外延迟或扰动、执行器能够施加计算量时，
 * 把这个 u[k] 代回对象方程，就得到 y[k+1] = r[k+1]。
 *
 * 本例先做手算检查，再比较理想死拍、限幅死拍和 PI，
 * 最后单独展示模型失配与测量噪声的影响。详细推导见同目录讲解.md。
 *
 * 本文件把算法和教学检查放在一起。main 中的 assert 会执行函数调用，
 * 验证时不能定义 NDEBUG；移植到工程时应显式调用并检查返回值。
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "../公共代码/PID核心.h"

/* 模型系数已经包含采样周期的影响，计算时不要再随意乘一次 dt。 */
typedef struct {
    double a;        /* 状态保留系数，无量纲：上一拍输出留下多少。 */
    double b;        /* 输入增益，单位为“被控输出单位/控制输入单位”。 */
    double minimum;  /* 执行器实际允许的控制输入下限，与 u 同单位。 */
    double maximum;  /* 执行器实际允许的控制输入上限，与 u 同单位。 */
} Model;

/* 区分“控制器想要多少”和“执行器实际能给多少”，便于观察限幅。 */
typedef struct {
    double requested; /* 由模型反解的理想控制请求，与 u 同单位。 */
    double applied;   /* 限幅后的实际控制输入；应把这个量送给对象。 */
    int limited;      /* 1 表示本次请求被截断，0 表示请求能够完整施加。 */
} Command;

/*
 * 计算一拍控制量，本函数不更新对象，也不保存控制历史。
 *
 * m：控制器使用的标称模型，不保证与真实对象完全相同。
 * measurement：当前测量 y[k]，与目标使用相同单位。
 * next_reference：下一拍目标 r[k+1]；本例为常数 0.8。
 * out：调用者提供的结果地址；成功后通过 *out 写回完整结果。
 *
 * 成功返回 1，失败返回 0。失败时不修改 out 中原有的内容。
 * 调用者必须检查返回值，不能把失败后保留的旧结果当作新命令。
 */
static int deadbeat(const Model *m, double measurement,
                    double next_reference, Command *out)
{
    /*
     * 先检查指针，再读取 m 的字段；|| 会短路，不会解引用空指针。
     * isfinite 拒绝 NaN 和无穷大。b 太小会造成巨大请求或除零风险。
     * 1e-12 是按本例归一化尺度选择的阈值，不是通用物理限值。
     * b 可以为负；只要方向正确且幅值足够，仍然可以反解。
     */
    if (!m || !out || !isfinite(m->a) || !isfinite(m->b) ||
        fabs(m->b) < 1e-12 ||
        !isfinite(m->minimum) || !isfinite(m->maximum) ||
        m->minimum >= m->maximum ||
        !isfinite(measurement) || !isfinite(next_reference))
        return 0;

    /* 第一步：求理想请求。m->a 表示读取模型指针所指向对象的 a。 */
    Command next;
    next.requested = (next_reference - m->a * measurement) / m->b;
    if (!isfinite(next.requested))
        return 0;

    /*
     * 第二步：限制到执行器范围。
     * 限幅后通常不再满足理想求解方程，因此“一拍到达”可能失效。
     * control_clamp 复用公共 PID 头文件中的限幅函数。
     */
    next.applied = control_clamp(next.requested, m->minimum, m->maximum);
    next.limited = (next.applied != next.requested);

    /* 第三步：确认计算有效后统一写回，避免失败时只更新部分字段。 */
    *out = next;
    return 1;
}

/* 一次完整仿真的汇总指标；本例 y、u 都使用归一化数值。 */
typedef struct {
    double iae;        /* 绝对误差积分近似值，单位为“被控输出单位*秒”。 */
    double final;      /* 最后一拍的对象输出，与 y 同单位。 */
    int first_near;    /* 首次进入目标 ±1% 的拍数；-1 表示始终未进入。 */
    int limited_steps; /* 死拍分支触发限幅的拍数；PI 分支未统计此项。 */
} Result;

/*
 * use_deadbeat 非零时使用死拍，否则使用公共 PID 的 PI 配置。
 * limit 指定两条路径相同形式的输入范围 [-limit, limit]。
 * real_a、real_b 是仿真真实对象系数；控制器始终使用标称 a=0.9、b=0.1。
 * 因此修改 real_a 可以单独观察“模型估计不准”带来的后果。
 */
static Result simulate(int use_deadbeat, double limit,
                       double real_a, double real_b)
{
    const Model nominal = {0.9, 0.1, -limit, limit};

    /*
     * PidConfig 字段顺序：Kp、Ki、Kd、微分滤波时间常数、下限、上限。
     * Kd=0，所以这里是 PI；公共核心还带有条件积分抗饱和。
     * 这组参数用于建立对照，不声称是针对该对象的最优整定结果。
     */
    const PidConfig cfg = {1.5, 2, 0, 0, -limit, limit};
    PidState pi = {0}; /* 每组仿真都从零积分和未初始化的状态开始。 */

    const double dt = 0.1, reference = 0.8; /* 采样周期 0.1 秒，恒定目标 0.8。 */
    double y = 0, area = 0;
    int first = -1, saturated = 0;

    /* 200 拍对应 20 秒。每一拍按“测量→控制→对象更新→统计”执行。 */
    for (int k = 0; k < 200; k++) {
        double u;

        /* 1. 使用当前 y[k] 计算本拍输入 u[k]。 */
        if (use_deadbeat) {
            Command command;
            assert(deadbeat(&nominal, y, reference, &command));
            u = command.applied;
            saturated += command.limited;
        } else {
            /* 参数 0 是前馈输入，本例没有额外前馈。 */
            assert(pid_update(&pi, &cfg, reference, y, 0, dt, &u));
        }
        assert(fabs(u) <= limit);

        /* 2. 真实对象从 y[k] 走到 y[k+1]，只接收实际施加量 u。 */
        y = real_a * y + real_b * u;
        assert(isfinite(y));

        /*
         * 3. 用更新后的 y[k+1] 统计误差。
         * area 累加 |目标-输出|*dt；它没有统计采样点之间的连续轨迹。
         * 所以理想死拍的 area 接近 0，不表示初始瞬间也没有误差。
         */
        area += fabs(reference - y) * dt;

        /*
         * 本例目标为正数 0.8，1% 误差带的半宽是 0.008。
         * k 从 0 开始，完成一次对象更新后已经是第 k+1 拍。
         * first<0 保证只记录首次进入，不能把它称作永不离开的整定时间。
         */
        if (first < 0 && fabs(reference - y) <= 0.01 * reference)
            first = k + 1;
    }

    Result result = {area, y, first, saturated};
    return result;
}

int main(void)
{
    /*
     * 检查一：手算一拍。
     * y=0.2、目标=0.8 时，u=(0.8-0.9*0.2)/0.1=6.2。
     * ±100 的范围足以施加 6.2，代回真实模型后下一拍应为 0.8。
     */
    Model m = {0.9, 0.1, -100, 100};
    Command c;
    assert(deadbeat(&m, 0.2, 0.8, &c));
    assert(fabs(c.requested - 6.2) < 1e-12 && !c.limited);
    assert(fabs(0.9 * 0.2 + 0.1 * c.applied - 0.8) < 1e-12);

    /* 负输入增益同样可解：b=-0.1 时需要 -6.2，不能擅自取绝对值。 */
    Model negative = {0.9, -0.1, -100, 100};
    assert(deadbeat(&negative, 0.2, 0.8, &c) &&
           fabs(c.applied + 6.2) < 1e-12);

    /*
     * 检查二：限幅与非法参数。
     * 从 y=0 到目标 0.8，理想请求是 8；执行器最多只能给 1。
     * 随后把 b 改成 0，验证求解失败且已有请求量、施加量保持不变。
     */
    m.maximum = 1;
    m.minimum = -1;
    assert(deadbeat(&m, 0, 0.8, &c) && c.limited && c.applied == 1);

    Command saved = c;
    m.b = 0;
    assert(!deadbeat(&m, 0, 0.8, &c) &&
           c.applied == saved.applied && c.requested == saved.requested);

    /*
     * 检查三：四组闭环。
     * ideal：输入范围足够大，用于验证一拍到达的数学性质。
     * limited 与 pi：同为 ±1，比较相同执行器能力下的结果。
     * mismatch：真实 a=0.85，仍给 ±100，单独观察模型失配。
     */
    Result ideal = simulate(1, 100, 0.9, 0.1);
    Result limited = simulate(1, 1, 0.9, 0.1);
    Result pi = simulate(0, 1, 0.9, 0.1);
    Result mismatch = simulate(1, 100, 0.85, 0.1);

    assert(ideal.first_near == 1 && ideal.iae < 1e-12);
    assert(limited.first_near > 1 && limited.limited_steps > 0);
    assert(limited.iae < pi.iae && fabs(limited.final - 0.8) < 1e-12);
    assert(fabs(pi.final - 0.8) < 1e-6);

    /*
     * 失配时 y_next=(0.85-0.9)*y+0.8=-0.05*y+0.8。
     * 稳态满足 1.05*y=0.8，因此最终值是 0.8/1.05，存在稳态误差。
     */
    assert(fabs(mismatch.final - 0.8 / 1.05) < 1e-12);

    /*
     * 检查四：测量噪声。
     * 恢复足够大的输出范围，仅把测量从 0.2 改成 0.21。
     * 控制量变化为 -(a/b)*0.01=-0.09，说明模型逆计算会放大测量误差。
     */
    m = (Model){0.9, 0.1, -100, 100};
    Command clean, noisy;
    assert(deadbeat(&m, 0.2, 0.8, &clean) &&
           deadbeat(&m, 0.21, 0.8, &noisy));
    assert(fabs((noisy.applied - clean.applied) + 0.09) < 1e-12);

    /* 打印固定模型下的实际指标；详细解释与适用范围见讲解.md。 */
    printf("deadbeat: first 1%% band step ideal=%d limited=%d PI=%d; limited steps=%d\n",
           ideal.first_near, limited.first_near, pi.first_near, limited.limited_steps);
    printf("deadbeat: 20s IAE limited=%.6f PI=%.6f; mismatched final=%.6f\n",
           limited.iae, pi.iae, mismatch.final);
    return 0;
}

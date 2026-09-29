/*
 * FOC的输入和运行管理辅助函数，不访问ADC、编码器寄存器或栅极硬件。
 * “门极禁止”必须由调用者真正落实到驱动硬件，不能只输出50%占空比。
 */
#ifndef FOC_INPUT_GUARD_H
#define FOC_INPUT_GUARD_H
#include "电机坐标变换.h"
#include "PID核心.h"

#define FOC_ZERO_SAMPLES 64
typedef struct {
    MotorAbc offset; /* 已换算成A的三路平均零偏，不是ADC计数。 */
    int count;      /* 成功采集次数；结构启动时用{0}。 */
    int ready;
} FocZero;

/*
 * 只有确认电流为零且门极禁止时采集。
 * 仅门极禁止并不能保证电机电流已经消失，zero_confirmed需由外部确认。
 * 不在有电流的运行过程中“自动校零”，避免把真实电流当成零偏。
 */
static inline int foc_zero_sample(FocZero *s, MotorAbc sample,
                                  int gates_disabled, int zero_confirmed)
{
    FocZero next;
    if (!s || s->ready || s->count < 0 || s->count >= FOC_ZERO_SAMPLES ||
        !gates_disabled || !zero_confirmed ||
        !isfinite(sample.a) || !isfinite(sample.b) || !isfinite(sample.c))
        return 0;
    next = *s;
    ++next.count;
    /* 增量均值，减少大数组存储；next仅在全部有限时写回。 */
    next.offset.a += (sample.a - next.offset.a) / next.count;
    next.offset.b += (sample.b - next.offset.b) / next.count;
    next.offset.c += (sample.c - next.offset.c) / next.count;
    if (!isfinite(next.offset.a) || !isfinite(next.offset.b) || !isfinite(next.offset.c))
        return 0;
    next.ready = next.count == FOC_ZERO_SAMPLES;
    *s = next;
    return 1;
}

static inline int foc_zero_correct(const FocZero *s, MotorAbc raw, MotorAbc *out)
{
    MotorAbc next;
    if (!s || !s->ready || !out)
        return 0;
    next.a = raw.a - s->offset.a;
    next.b = raw.b - s->offset.b;
    next.c = raw.c - s->offset.c;
    if (!isfinite(next.a) || !isfinite(next.b) || !isfinite(next.c))
        return 0;
    *out = next;
    return 1;
}

/* 零位偏置按电角度定义；direction必须±1，速度反馈也要采用相同方向。 */
static inline int foc_electrical_angle(double mechanical, int pairs,
                                       int direction, double offset, double *out)
{
    if (!out || !isfinite(mechanical) || !isfinite(offset) ||
        pairs < 1 || pairs > 64 || (direction != 1 && direction != -1))
        return 0;
    double angle = direction * pairs * mechanical + offset;
    if (!isfinite(angle))
        return 0;
    angle = fmod(angle, 2.0 * MOTOR_PI);
    *out = angle < 0.0 ? angle + 2.0 * MOTOR_PI : angle;
    return 1;
}

/* 机械速度参考斜坡：rate是rad/s^2，dt是速度环周期秒数。 */
static inline int foc_ramp(double current, double target, double rate,
                           double dt, double *out)
{
    if (!out || !isfinite(current) || !isfinite(target) ||
        !isfinite(rate) || rate < 0.0 || !isfinite(dt) || dt <= 0.0)
        return 0;
    double change = target - current, step = rate * dt;
    if (!isfinite(change) || !isfinite(step))
        return 0;
    double next = current + control_clamp(change, -step, step);
    if (!isfinite(next))
        return 0;
    *out = next;
    return 1;
}

/*
 * d轴优先的电流参考圆：先限制id，再按剩余余量限制iq。
 * 这是参考约束，不是实际电流保护；实际电流超限仍由guard处理。
 */
static inline int foc_current_circle(double id, double iq, double limit,
                                     MotorDq0 *out)
{
    if (!out || !isfinite(id) || !isfinite(iq) ||
        !isfinite(limit) || limit <= 0.0)
        return 0;
    MotorDq0 next = {control_clamp(id, -limit, limit), 0.0, 0.0};
    double ratio = next.d / limit;
    double q_limit = limit * sqrt(fmax(0.0, 1.0 - ratio * ratio));
    next.q = control_clamp(iq, -q_limit, q_limit);
    *out = next;
    return 1;
}

enum {
    FOC_BAD_INPUT = 1u,
    FOC_UNDERVOLTAGE = 2u,
    FOC_OVERVOLTAGE = 4u,
    FOC_OVERCURRENT = 8u,
    FOC_OVERSPEED = 16u
};
typedef struct {
    double bus_min, bus_max; /* 母线允许范围[V]。 */
    double phase_current;   /* 任意一相电流绝对值跳闸阈值[A]。 */
    double speed;           /* 机械角速度绝对值跳闸阈值[rad/s]。 */
} FocGuardLimits;
typedef struct {
    unsigned faults; /* 位掩码锁存，可同时有多个故障。 */
    int gate_enable; /* 1允许驱动，0必须禁止门极；不是某个PWM占空比。 */
} FocGuard;

/*
 * 每个电流周期先检查，再决定是否允许控制。
 * 故障恢复不自动清除；仅停机请求且当前检查健康时接受reset。
 * 未完成零偏校准只禁止启动，不把等待校准本身作为硬件故障。
 */
static inline int foc_guard_update(FocGuard *s, FocGuardLimits cfg,
                                   int enable, int reset, MotorAbc current,
                                   double bus, double speed, int angle_valid,
                                   int calibrated)
{
    unsigned active = 0u;
    if (!s)
        return 0;
    if (!isfinite(cfg.bus_min) || !isfinite(cfg.bus_max) ||
        cfg.bus_min <= 0.0 || cfg.bus_min >= cfg.bus_max ||
        !isfinite(cfg.phase_current) || cfg.phase_current <= 0.0 ||
        !isfinite(cfg.speed) || cfg.speed <= 0.0 ||
        !isfinite(bus) || !isfinite(speed) ||
        !isfinite(current.a) || !isfinite(current.b) || !isfinite(current.c) ||
        !angle_valid) {
        active = FOC_BAD_INPUT;
    } else {
        if (bus < cfg.bus_min) active |= FOC_UNDERVOLTAGE;
        if (bus > cfg.bus_max) active |= FOC_OVERVOLTAGE;
        if (fabs(current.a) > cfg.phase_current ||
            fabs(current.b) > cfg.phase_current ||
            fabs(current.c) > cfg.phase_current) active |= FOC_OVERCURRENT;
        if (fabs(speed) > cfg.speed) active |= FOC_OVERSPEED;
    }
    s->faults |= active;
    if (reset && !enable && active == 0u)
        s->faults = 0u;
    s->gate_enable = enable && calibrated && s->faults == 0u;
    return s->gate_enable;
}

/*
 * 只支持PI（Kd=0）。内环电压不足且积分想加重缺额时，冻结外环积分。
 * deficit=上一电流拍vq_raw-vq_applied[V]；使用上一拍反馈，存在一拍延迟。
 * error与deficit符号相反时，仍允许积分帮助恢复。
 */
static inline int foc_speed_pi(PidState *s, const PidConfig *cfg,
                               double reference, double measured, double dt,
                               int limited, double deficit, double *out)
{
    PidState next;
    double command;
    if (!s || !out || !pid_config_valid(cfg) || cfg->kd != 0.0 ||
        !isfinite(deficit))
        return 0;
    next = *s;
    if (!pid_update(&next, cfg, reference, measured, 0.0, dt, &command))
        return 0;
    double error = reference - measured;
    if (limited && ((error > 0.0 && deficit > 0.0) ||
                    (error < 0.0 && deficit < 0.0))) {
        next.integral = cfg->ki == 0.0 ? 0.0 : s->integral;
        double raw = cfg->kp * error + next.integral;
        if (!isfinite(raw))
            return 0;
        command = control_clamp(raw, cfg->minimum, cfg->maximum);
    }
    *s = next;
    *out = command;
    return 1;
}
#endif

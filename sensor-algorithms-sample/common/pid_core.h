/*
 * 位置式 PID 公共核心：
 * u = feedforward + Kp*error + integral - Kd*filtered_measurement_rate。
 * 输出限幅、条件积分抗饱和、测量微分低通；不依赖特定硬件。
 */
#ifndef SENSOR_PID_CORE_H
#define SENSOR_PID_CORE_H
#include <math.h>
#include <stddef.h>

typedef struct
{
    double kp, ki, kd;        /* 比例、积分、微分系数；与 dt 秒配套。 */
    double derivative_tau;  /* 测量微分低通时间常数，秒；0 表示不过滤。 */
    double minimum, maximum;/* 实际执行器能够接受的输出范围。 */
} PidConfig;
typedef struct
{
    double integral;        /* 已乘 Ki 的积分贡献，单位与控制输出相同。 */
    double previous_measurement;
    double derivative;      /* 滤波后的测量变化率。 */
    int initialized;
} PidState;

static inline double control_clamp(double value, double minimum, double maximum)
{
    return fmax(minimum,fmin(maximum,value));
}
static inline int pid_config_valid(const PidConfig *cfg)
{
    return cfg!=NULL && isfinite(cfg->kp) && cfg->kp>=0 &&
        isfinite(cfg->ki) && cfg->ki>=0 && isfinite(cfg->kd) && cfg->kd>=0 &&
        isfinite(cfg->derivative_tau) && cfg->derivative_tau>=0 &&
        isfinite(cfg->minimum) && isfinite(cfg->maximum) && cfg->minimum<cfg->maximum;
}
/*
 * state 初始用 {0}。output 不得指向 state 或 cfg 的成员。
 * 输入要求正作用对象：控制输出增加时，被控测量通常也增加。
 * 返回 0 时不改状态和输出；返回 1 表示本次计算成功。
 */
static inline int pid_update(PidState *state, const PidConfig *cfg,
                              double reference, double measurement,
                              double feedforward, double dt, double *output)
{
    if (state==NULL || output==NULL || !pid_config_valid(cfg) ||
        !isfinite(reference) || !isfinite(measurement) || !isfinite(feedforward) ||
        !isfinite(dt) || dt<=0 || !isfinite(state->integral) ||
        !isfinite(state->derivative) || !isfinite(state->previous_measurement))
        return 0;
    PidState next=*state;
    const double error=reference-measurement;
    /*
     * 对测量微分而不是对误差微分：设定值突然改变不会直接触发 D 冲击。
     * 首次无历史测量，导数取 0，不把默认的旧零值当真实历史。
     */
    double measured_rate=0;
    if (state->initialized)
        measured_rate=(measurement-state->previous_measurement)/dt;
    const double alpha=dt/(cfg->derivative_tau+dt);
    next.derivative=state->initialized
        ? state->derivative+alpha*(measured_rate-state->derivative) : 0;
    if (!isfinite(error) || !isfinite(measured_rate) || !isfinite(next.derivative))
        return 0;

    /* 先尝试积分，再判断这一步是否会让饱和更严重。 */
    const double previous_integral=cfg->ki==0 ? 0 : state->integral;
    const double increment=cfg->ki*error*dt;
    next.integral=previous_integral+increment;
    const double proportional=cfg->kp*error;
    const double differential=-cfg->kd*next.derivative;
    double raw=feedforward+proportional+next.integral+differential;
    if (!isfinite(next.integral) || !isfinite(raw)) return 0;

    /*
     * 上限饱和且积分还想增加，或者下限饱和且积分还想减少：
     * 拒绝这次积分。反方向积分仍允许，帮助退出饱和。
     * feedforward 一并计入 raw，避免忽略前馈已经占用的执行器余量。
     */
    if ((raw>cfg->maximum && increment>0) ||
        (raw<cfg->minimum && increment<0))
    {
        next.integral=previous_integral;
        raw=feedforward+proportional+next.integral+differential;
    }
    if (!isfinite(raw)) return 0;
    const double command=control_clamp(raw,cfg->minimum,cfg->maximum);
    next.previous_measurement=measurement;
    next.initialized=1;
    *state=next;
    *output=command;
    return 1;
}
#endif

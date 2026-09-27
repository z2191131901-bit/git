/*
 * 二阶 IIR 公共核心：系数已除以 a0。
 * y[n]=b0*x[n]+b1*x[n-1]+b2*x[n-2]-a1*y[n-1]-a2*y[n-2]。
 * 使用转置直接 II 型，仅需保存两个延迟状态 z1/z2。
 */
#ifndef SENSOR_BIQUAD_H
#define SENSOR_BIQUAD_H
#include "math3d.h"
typedef struct
{
    double b0,b1,b2,a1,a2;
    double z1,z2;
} Biquad;

static inline int biquad_design(Biquad *out, double fs, double frequency,
                                double quality, int notch)
{
    if (out==NULL || !isfinite(fs) || fs<=0 || !isfinite(frequency) ||
        frequency<=0 || frequency>=fs/2 || !isfinite(quality) || quality<=0)
        return 0;
    const double omega=2*SENSOR_PI*(frequency/fs);
    const double cosine=cos(omega);
    const double alpha=sin(omega)/(2*quality);
    const double a0=1+alpha;
    Biquad next={0};
    /* 低通分子为 (1-cos w)/2*[1,2,1]；陷波分子在目标频率产生零点。 */
    next.b0=(notch?1:(1-cosine)/2)/a0;
    next.b1=(notch?-2*cosine:1-cosine)/a0;
    next.b2=next.b0;
    next.a1=-2*cosine/a0;
    next.a2=(1-alpha)/a0;
    if (!isfinite(a0) || !isfinite(next.a1) || !isfinite(next.a2) ||
        1+next.a1+next.a2<=0 || 1-next.a1+next.a2<=0 || fabs(next.a2)>=1)
        return 0; /* 浮点精度下也必须满足二阶稳定性条件。 */
    *out=next; /* 重设计会清零状态，不能每个样本都重新设计。 */
    return 1;
}
static inline int biquad_update(Biquad *state, double input, double *output)
{
    if (state==NULL || output==NULL || !isfinite(input)) return 0;
    const double y=state->b0*input+state->z1;
    const double z1=state->b1*input-state->a1*y+state->z2;
    const double z2=state->b2*input-state->a2*y;
    if (!isfinite(y) || !isfinite(z1) || !isfinite(z2)) return 0;
    state->z1=z1; state->z2=z2; *output=y;
    return 1;
}
#endif

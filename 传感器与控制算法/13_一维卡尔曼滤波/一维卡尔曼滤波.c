/*
 * 一维卡尔曼：随机游走模型 x(k)=x(k-1)+过程噪声，z=x+测量噪声。
 * P/Q/R 都是方差（被测单位的平方），Q 是每步方差，不是标准差。
 */
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

typedef struct { double value, variance; } ScalarKalman;
static int kalman_update(ScalarKalman *state, double measurement, int valid,
                         double process_variance, double measurement_variance)
{
    if (state==NULL || !isfinite(state->value) || !isfinite(state->variance) ||
        state->variance<0 || !isfinite(process_variance) || process_variance<0 ||
        !isfinite(measurement_variance) || measurement_variance<=0 ||
        (valid && !isfinite(measurement))) return 0;
    ScalarKalman next=*state;
    /* ① 没有已知运动量，预测值不变，但预测方差增加。 */
    next.variance+=process_variance;
    if (!isfinite(next.variance)) return 0;
    if (valid)
    {
        /* ② 增益越大，测量对本次结果的影响越大。 */
        const double innovation_variance=next.variance+measurement_variance;
        if (!isfinite(innovation_variance)) return 0;
        const double gain=next.variance/innovation_variance;
        next.value+=gain*(measurement-next.value);
        /* ③ Joseph 形式标量协方差更新，保持非负的数值表达。 */
        next.variance=(1-gain)*(1-gain)*next.variance
                     +gain*gain*measurement_variance;
    }
    if (!isfinite(next.value) || !isfinite(next.variance)) return 0;
    *state=next;
    return 1;
}
int main(void)
{
    ScalarKalman state={10,4};
    int ok=kalman_update(&state,14,1,0,4);
    assert(ok && fabs(state.value-12)<1e-12 && fabs(state.variance-2)<1e-12);
    printf("estimate=%.3f variance=%.3f\n",state.value,state.variance);

    /* 测量丢失时只增加不确定性，不拿 NaN 当观测。 */
    ok=kalman_update(&state,NAN,0,1,4);
    assert(ok && state.value==12 && state.variance==3);
    ScalarKalman saved=state;
    ok=kalman_update(&state,14,1,1,0);
    assert(!ok && state.variance==saved.variance);
    for (int i=0;i<100;++i)
    {
        ok=kalman_update(&state,14,1,0.1,1);
        assert(ok && state.variance>=0);
    }
    assert(fabs(state.value-14)<1e-10);
    puts("kalman_scalar: PASS");
    return 0;
}

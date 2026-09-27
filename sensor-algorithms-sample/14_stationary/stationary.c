/*
 * 滑动窗口静止检测：32 个样本全部满足低角速度、接近 g，
 * 且加速度三轴方差之和足够小，才判为静止。
 * 单位：加速度 m/s^2，角速度 rad/s。不能检测匀速直线运动。
 */
#include "../common/math3d.h"
#include <assert.h>
#include <stdio.h>
#define STILL_WINDOW 32

typedef struct
{
    Vec3 acceleration[STILL_WINDOW];
    double gyro_norm[STILL_WINDOW];
    size_t next, count;
} Stationary;

static int stationary_update(Stationary *state, Vec3 accel, Vec3 gyro, int *is_still)
{
    if (state==NULL || is_still==NULL || !vec_finite(accel) || !vec_finite(gyro))
        return 0;
    const double angular_speed=vec_norm(gyro);
    if (!isfinite(angular_speed) || !isfinite(vec_norm(accel))) return 0;
    /* 先操作副本，溢出或异常时保留调用前状态。 */
    Stationary next=*state;
    next.acceleration[next.next]=accel;
    next.gyro_norm[next.next]=angular_speed;
    next.next=(next.next+1)%STILL_WINDOW;
    if (next.count<STILL_WINDOW) next.count++;
    if (next.count<STILL_WINDOW)
    {
        *state=next; *is_still=0; /* 尚未填满窗口，不能提前确认静止。 */
        return 1;
    }

    Vec3 mean={0,0,0};
    double m2=0;
    int within_limits=1;
    for (size_t i=0;i<STILL_WINDOW;++i)
    {
        Vec3 a=next.acceleration[i];
        if (fabs(vec_norm(a)-9.80665)>0.3 || next.gyro_norm[i]>0.03)
            within_limits=0;
        /*
         * Welford 在线方差：先计算旧均值差，再更新均值，
         * 最后用旧差乘新差，避免“大数平方相减”的精度损失。
         */
        Vec3 delta=vec_add(a,vec_scale(mean,-1));
        mean=vec_add(mean,vec_scale(delta,1.0/(double)(i+1)));
        Vec3 after=vec_add(a,vec_scale(mean,-1));
        m2+=delta.x*after.x+delta.y*after.y+delta.z*after.z;
    }
    if (!vec_finite(mean) || !isfinite(m2)) return 0;
    const double variance=m2/(STILL_WINDOW-1);
    *is_still=within_limits && variance<0.01;
    *state=next;
    return 1;
}
int main(void)
{
    Stationary state={0};
    Vec3 up={0,0,9.80665}, zero={0,0,0};
    int still=0,ok;
    for (int i=0;i<STILL_WINDOW;++i)
    {
        ok=stationary_update(&state,up,zero,&still);
        assert(ok && still==(i==STILL_WINDOW-1));
    }
    Vec3 turn={0,0,0.2};
    ok=stationary_update(&state,up,turn,&still);
    assert(ok && !still);
    for (int i=0;i<STILL_WINDOW;++i)
    {
        ok=stationary_update(&state,up,zero,&still);
        assert(ok);
    }
    assert(still);
    ok=stationary_update(&state,zero,zero,&still);
    assert(ok && !still); /* 自由落体/零读数不是静止。 */
    size_t saved=state.next;
    ok=stationary_update(&state,(Vec3){NAN,0,0},zero,&still);
    assert(!ok && state.next==saved);

    state=(Stationary){0};
    for (int i=0;i<STILL_WINDOW;++i)
    {
        /* 模长接近 g，但横向正负振动导致方差偏大。 */
        ok=stationary_update(&state,(Vec3){i%2?0.2:-0.2,0,9.80665},zero,&still);
        assert(ok);
    }
    assert(!still);
    puts("stationary: PASS");
    return 0;
}

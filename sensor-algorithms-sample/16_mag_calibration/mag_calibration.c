/*
 * 磁力计 min/max 校准：硬铁偏移 + 轴对齐的独立比例校正。
 * 不是完整椭球拟合，不能修正任意非对角软铁耦合。
 * 输入同一磁场单位；输出归一到参考球半径 1（无量纲）。
 */
#include "../common/math3d.h"
#include <assert.h>
#include <stdio.h>

typedef struct { Vec3 bias, radius; } MagCalibration;
static int mag_fit(const Vec3 *samples, size_t count, MagCalibration *out)
{
    if (samples==NULL || out==NULL || count<6) return 0;
    if (!vec_finite(samples[0])) return 0;
    Vec3 low=samples[0],high=samples[0];
    for (size_t i=1;i<count;++i)
    {
        Vec3 v=samples[i];
        if (!vec_finite(v)) return 0;
        low.x=fmin(low.x,v.x); low.y=fmin(low.y,v.y); low.z=fmin(low.z,v.z);
        high.x=fmax(high.x,v.x); high.y=fmax(high.y,v.y); high.z=fmax(high.z,v.z);
    }
    MagCalibration next;
    next.bias=vec_add(vec_scale(high,0.5),vec_scale(low,0.5));
    next.radius=vec_add(vec_scale(high,0.5),vec_scale(low,-0.5));
    if (!vec_finite(next.radius) || next.radius.x<1e-9 ||
        next.radius.y<1e-9 || next.radius.z<1e-9) return 0;
    *out=next;
    return 1;
}
static int mag_correct(const MagCalibration *cal, Vec3 raw, Vec3 *out)
{
    if (cal==NULL || out==NULL || !vec_finite(raw) || !vec_finite(cal->bias) ||
        !vec_finite(cal->radius) || cal->radius.x<=0 || cal->radius.y<=0 ||
        cal->radius.z<=0) return 0;
    Vec3 next={(raw.x-cal->bias.x)/cal->radius.x,
               (raw.y-cal->bias.y)/cal->radius.y,
               (raw.z-cal->bias.z)/cal->radius.z};
    if (!vec_finite(next)) return 0;
    *out=next;
    return 1;
}
int main(void)
{
    /* 球经过偏移 [10,-5,3] 与独立伸缩 [2,3,4]，得到轴对齐椭球。 */
    const Vec3 samples[]={{12,-5,3},{8,-5,3},{10,-2,3},
                          {10,-8,3},{10,-5,7},{10,-5,-1}};
    MagCalibration cal;
    int ok=mag_fit(samples,6,&cal);
    assert(ok && cal.bias.x==10 && cal.radius.z==4);
    Vec3 diagonal={10+2/sqrt(3.0),-5+3/sqrt(3.0),3+4/sqrt(3.0)}, corrected;
    ok=mag_correct(&cal,diagonal,&corrected);
    assert(ok && fabs(vec_norm(corrected)-1)<1e-12);
    const Vec3 bad[6]={{0,0,0},{1,0,0},{2,0,0},{3,0,0},{4,0,0},{5,0,0}};
    ok=mag_fit(bad,6,&cal);
    assert(!ok); /* 只沿一条轴采集，无法标定三维范围。 */
    puts("mag_calibration: PASS");
    return 0;
}

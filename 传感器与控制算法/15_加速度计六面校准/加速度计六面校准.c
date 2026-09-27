/*
 * 简化六面校准：各轴独立的比例因子+零偏，不估计轴间耦合。
 * positive[i] 是第 i 轴朝上时该轴多次读数的平均，negative 为朝下。
 * 输入输出 m/s^2；每个面必须真正静止且对应轴准确对齐竖直。
 */
#include "../公共代码/三维数学.h"
#include <assert.h>
#include <stdio.h>

typedef struct { double bias[3], scale[3]; } AccelCalibration;
static int accel_calibrate(const double positive[3], const double negative[3],
                            double gravity, AccelCalibration *output)
{
    if (positive==NULL || negative==NULL || output==NULL ||
        !isfinite(gravity) || gravity<=0) return 0;
    AccelCalibration next;
    for (size_t axis=0;axis<3;++axis)
    {
        if (!isfinite(positive[axis]) || !isfinite(negative[axis])) return 0;
        const double span=positive[axis]-negative[axis];
        if (!isfinite(span) || span<=1e-9) return 0;
        /* 测量模型 raw=scale*true+bias；校正时使用 (raw-bias)/scale。 */
        next.bias[axis]=positive[axis]/2+negative[axis]/2;
        next.scale[axis]=span/(2*gravity);
        if (!isfinite(next.scale[axis]) || next.scale[axis]<=0) return 0;
    }
    *output=next;
    return 1;
}
static int accel_correct(const AccelCalibration *cal, const double raw[3], double out[3])
{
    if (cal==NULL || raw==NULL || out==NULL) return 0;
    double next[3];
    for (size_t i=0;i<3;++i)
    {
        if (!isfinite(raw[i]) || !isfinite(cal->bias[i]) ||
            !isfinite(cal->scale[i]) || cal->scale[i]<=0) return 0;
        next[i]=(raw[i]-cal->bias[i])/cal->scale[i];
        if (!isfinite(next[i])) return 0;
    }
    for (size_t i=0;i<3;++i) out[i]=next[i];
    return 1;
}
int main(void)
{
    const double g=9.80665, bias[3]={0.2,-0.1,0.3}, scale[3]={1.1,0.9,1.05};
    double positive[3],negative[3],raw[3],corrected[3];
    for (size_t i=0;i<3;++i)
    {
        positive[i]=bias[i]+scale[i]*g;
        negative[i]=bias[i]-scale[i]*g;
        raw[i]=bias[i]+scale[i]*(double)(i+1);
    }
    AccelCalibration cal;
    int ok=accel_calibrate(positive,negative,g,&cal);
    assert(ok);
    ok=accel_correct(&cal,raw,corrected);
    assert(ok);
    for (size_t i=0;i<3;++i)
        assert(fabs(corrected[i]-(double)(i+1))<1e-12);
    AccelCalibration saved=cal;
    positive[0]=negative[0];
    ok=accel_calibrate(positive,negative,g,&cal);
    assert(!ok && cal.scale[0]==saved.scale[0]);
    puts("accel_calibration: PASS");
    return 0;
}

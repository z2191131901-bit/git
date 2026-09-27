/* 二阶 Butterworth 低通：Q=1/sqrt(2)，截止频率处幅值约为 0.707。 */
#include "../公共代码/二阶滤波器.h"
#include <assert.h>
#include <stdio.h>

/* 生成正弦波并计算稳态输出 RMS，以验证真实频率响应。 */
static double response(double signal_frequency)
{
    Biquad filter;
    int ok=biquad_design(&filter,1000,50,1/sqrt(2.0),0);
    assert(ok);
    double energy=0,output=0;
    for (int i=0;i<4000;++i)
    {
        double input=sin(2*SENSOR_PI*signal_frequency*i/1000);
        ok=biquad_update(&filter,input,&output);
        assert(ok);
        if (i>=2000) energy+=output*output; /* 丢弃启动瞬态。 */
    }
    return sqrt(energy/2000)*sqrt(2.0); /* 输入正弦幅值为 1。 */
}
int main(void)
{
    const double cutoff_gain=response(50);
    const double stop_gain=response(400);
    assert(fabs(cutoff_gain-1/sqrt(2.0))<1e-8);
    assert(stop_gain<0.01);
    Biquad filter;
    int ok=biquad_design(&filter,1000,50,1/sqrt(2.0),0);
    assert(ok);
    double output=0;
    for (int i=0;i<1000;++i)
    {
        ok=biquad_update(&filter,1,&output);
        assert(ok);
    }
    assert(fabs(output-1)<1e-10);
    const double saved=filter.z1;
    ok=biquad_update(&filter,NAN,&output);
    assert(!ok && filter.z1==saved);
    ok=biquad_design(&filter,1000,500,1/sqrt(2.0),0);
    assert(!ok);
    printf("gain@50Hz=%.6f, gain@400Hz=%.6f\n",cutoff_gain,stop_gain);
    puts("butterworth: PASS");
    return 0;
}

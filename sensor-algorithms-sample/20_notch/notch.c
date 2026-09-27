/* 二阶陷波：抑制 50Hz 附近的窄带干扰，采样率 1000Hz，Q=10。 */
#include "../common/biquad.h"
#include <assert.h>
#include <stdio.h>

static double response(double frequency)
{
    Biquad filter;
    int ok=biquad_design(&filter,1000,50,10,1);
    assert(ok);
    double energy=0,output=0;
    for (int i=0;i<6000;++i)
    {
        ok=biquad_update(&filter,sin(2*SENSOR_PI*frequency*i/1000),&output);
        assert(ok);
        if (i>=4000) energy+=output*output;
    }
    return sqrt(energy/2000)*sqrt(2.0);
}
int main(void)
{
    const double rejected=response(50),passed=response(10);
    assert(rejected<1e-8 && passed>0.99 && passed<=1.000001);
    Biquad filter;
    int ok=biquad_design(&filter,1000,50,0,1);
    assert(!ok); /* Q 必须为正数。 */
    ok=biquad_design(&filter,1000,50,10,1);
    assert(ok);
    double output=0;
    for (int i=0;i<4000;++i)
    {
        ok=biquad_update(&filter,1,&output);
        assert(ok);
    }
    assert(fabs(output-1)<1e-10); /* 直流分量应通过。 */
    printf("gain@50Hz=%.9f, gain@10Hz=%.6f\n",rejected,passed);
    puts("notch: PASS");
    return 0;
}

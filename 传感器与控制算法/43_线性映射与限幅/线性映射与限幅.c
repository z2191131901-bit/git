/* 输入区间必须递增，输出区间可递增也可递减；区间外饱和，不外推。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
static int map_clamped(double x,double xmin,double xmax,
                       double ymin,double ymax,double *output)
{
    if(!output || !isfinite(x) || !isfinite(xmin) || !isfinite(xmax) ||
       !isfinite(ymin) || !isfinite(ymax) || xmin>=xmax) return 0;
    double span=xmax-xmin;
    if(!isfinite(span)) return 0;
    /* 先限幅，可避免远离区间的输入使 x-xmin 溢出。 */
    double clipped=fmax(xmin,fmin(xmax,x));
    double ratio=(clipped-xmin)/span;
    double result=(1-ratio)*ymin+ratio*ymax;
    if(!isfinite(result)) return 0;
    *output=result;
    return 1;
}
int main(void)
{
    double y=0;
    assert(map_clamped(2047.5,0,4095,0,3.3,&y) && fabs(y-1.65)<1e-12);
    assert(map_clamped(-100,0,4095,0,3.3,&y) && y==0);
    assert(map_clamped(5000,0,4095,0,3.3,&y) && y==3.3);
    assert(map_clamped(.25,0,1,100,0,&y) && y==75);
    assert(!map_clamped(0,1,1,0,1,&y) && y==75);
    puts("linear map: PASS, ADC midpoint=1.65 V");
    return 0;
}

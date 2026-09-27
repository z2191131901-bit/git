/* 峰值：严格三点局部极大值，延迟一个采样点确认，不检测平顶峰。
 * 过零：只统计负到非负的上升过零，线性插值估计过零时间。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
typedef struct {
    double older,previous,previous_time,last_cross;
    unsigned samples;
    int have_cross;
} Detector;
typedef struct { int peak,cross,frequency_valid; double peak_value,frequency; } Event;
static int update(Detector *d,double x,double time,double threshold,Event *event)
{
    if(!d || !event || !isfinite(x) || !isfinite(time) || !isfinite(threshold) ||
       (d->samples>0 && time<=d->previous_time)) return 0;
    Event e={0,0,0,0,0};
    if(d->samples>=2 && d->previous>d->older && d->previous>x &&
       d->previous>=threshold) { e.peak=1; e.peak_value=d->previous; }
    if(d->samples>0 && d->previous<0 && x>=0) {
        double fraction=-d->previous/(x-d->previous);
        double crossing=d->previous_time+fraction*(time-d->previous_time);
        if(!isfinite(crossing)) return 0;
        e.cross=1;
        if(d->have_cross) {
            double period=crossing-d->last_cross;
            if(period<=0 || !isfinite(period)) return 0;
            e.frequency=1/period;
            if(!isfinite(e.frequency)) return 0;
            e.frequency_valid=1;
        }
        d->last_cross=crossing; d->have_cross=1;
    }
    d->older=d->previous; d->previous=x; d->previous_time=time;
    if(d->samples<2) d->samples++;
    *event=e;
    return 1;
}
int main(void)
{
    Detector d={0};
    Event e={0};
    int peaks=0,crosses=0,frequencies=0;
    /* 2 Hz，100 Hz采样，避开正弦峰值恰好落在相邻等高点的情况。 */
    for(int i=0;i<=200;i++) {
        double t=i*.01;
        assert(update(&d,sin(2*3.14159265358979323846*2*t-.3),t,.8,&e));
        peaks+=e.peak; crosses+=e.cross;
        if(e.frequency_valid) { frequencies++; assert(fabs(e.frequency-2)<1e-10); }
    }
    assert(peaks==4 && crosses==4 && frequencies==3);
    assert(!update(&d,0,2,.8,&e));
    Detector plateau={0};
    const double flat[]={0,1,1,0};
    for(int i=0;i<4;i++) { assert(update(&plateau,flat[i],i,.5,&e)); assert(!e.peak); }
    puts("peak/zero crossing: PASS, four peaks, frequency=2 Hz");
    return 0;
}

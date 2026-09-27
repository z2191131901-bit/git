/* 仅模型上的继电反馈整定：对称继电输出 -> 周期/幅值 -> 近似Ku/Tu -> ZN PID。
 * 不是可直接连接硬件的自动整定器；没有传感器驱动和现场停机链。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "../公共代码/PID核心.h"
typedef struct { double x1,x2,y; } Plant;
/* G(s)=1/(s+1)^3，三阶同极点系统在零阶保持输入下精确传播。 */
static void plant(Plant *p,double u,double dt)
{
    double e=exp(-dt),x1=p->x1,x2=p->x2,y=p->y;
    p->x1=e*x1+(1-e)*u;
    p->x2=e*(x2+dt*x1)+(1-e*(1+dt))*u;
    p->y=e*(y+dt*x2+dt*dt*x1/2)+(1-e*(1+dt+dt*dt/2))*u;
}
typedef struct { double amplitude,period,ku,kp,ki,kd; int cycles; } Tuning;
static int tune(double height,double dt,Tuning *out)
{
    if(!out || !isfinite(height) || height<=0 || !isfinite(dt) || dt<.0001 || dt>.01)
        return 0;
    Plant p={0};
    double last_cross=0,previous_period=0,previous_amplitude=0,period_sum=0,amplitude_sum=0;
    double low=0,high=0;
    int have_cross=0,cycles=0;
    /* 丢弃前60秒暂态，最多仿真120秒；收集10个稳定周期后结束。 */
    for(int k=0;k<(int)(120/dt);k++) {
        double old=p.y,u=p.y<=0 ? height : -height;
        plant(&p,u,dt);
        if(!isfinite(p.y) || fabs(p.y)>2*height) return 0;
        low=fmin(low,p.y); high=fmax(high,p.y);
        if(old<=0 && p.y>0) {
            double crossing=k*dt+dt*(-old)/(p.y-old);
            if(have_cross && crossing>60) {
                double period=crossing-last_cross,amplitude=(high-low)/2;
                if(period<=0 || amplitude<1e-8) return 0;
                if(previous_period>0 && fabs(period-previous_period)>.05*previous_period)
                    return 0; /* 未形成近似稳定周期就不输出参数。 */
                if(previous_amplitude>0 && fabs(amplitude-previous_amplitude)>.05*previous_amplitude) return 0;
                period_sum+=period; amplitude_sum+=amplitude; previous_period=period; previous_amplitude=amplitude;
                cycles++;
                if(cycles==10) {
                    Tuning result={0};
                    result.period=period_sum/cycles;
                    result.amplitude=amplitude_sum/cycles;
                    result.ku=4*height/(3.14159265358979323846*result.amplitude);
                    result.kp=.6*result.ku;
                    double ti=.5*result.period,td=.125*result.period;
                    result.ki=result.kp/ti; result.kd=result.kp*td;
                    result.cycles=cycles; *out=result; return 1;
                }
            }
            last_cross=crossing; have_cross=1; low=p.y; high=p.y;
        }
    }
    return 0;
}
static double closed_loop(PidConfig cfg,double *final)
{
    Plant p={0}; PidState s={0};
    double area=0;
    for(int k=0;k<20000;k++) {
        double u;
        assert(pid_update(&s,&cfg,1,p.y,0,.002,&u));
        plant(&p,u,.002);
        assert(isfinite(p.y) && fabs(p.y)<5);
        area+=fabs(1-p.y)*.002;
    }
    *final=p.y; return area;
}
int main(void)
{
    Tuning t;
    assert(tune(.5,.002,&t));
    /* 该对象解析临界值Ku=8、Tu=2*pi/sqrt(3)。继电描述函数仅近似。 */
    double exact_period=2*3.14159265358979323846/sqrt(3);
    assert(fabs(t.ku-8)/8<.2 && fabs(t.period-exact_period)/exact_period<.2);
    assert(t.cycles==10 && t.kp>0 && t.ki>0 && t.kd>0);
    Tuning saved=t;
    assert(!tune(0,.002,&t) && t.ku==saved.ku);
    double baseline_end,tuned_end;
    double baseline_area=closed_loop((PidConfig){1,.2,0,.05,-5,5},&baseline_end);
    double tuned_area=closed_loop((PidConfig){t.kp,t.ki,t.kd,.05,-5,5},&tuned_end);
    assert(fabs(tuned_end-1)<.01 && tuned_area<baseline_area);
    printf("relay tuning: cycles=%d Ku=%.6f Tu=%.6f Kp=%.6f Ki=%.6f Kd=%.6f\n",
           t.cycles,t.ku,t.period,t.kp,t.ki,t.kd);
    printf("relay tuning validation: 40s IAE baseline=%.6f tuned=%.6f final=%.6f\n",
           baseline_area,tuned_area,tuned_end);
    return 0;
}

/* 一阶惯性+整数拍输入延迟的Smith预估器。
 * 反馈量=无延迟模型输出+(实际输出-带延迟模型输出)。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "../公共代码/PID核心.h"
enum { CAPACITY=200 };
typedef struct { double data[CAPACITY]; int index,length; } Delay;
static int delay_init(Delay *d,int length)
{
    if(!d || length<1 || length>CAPACITY) return 0;
    *d=(Delay){{0},0,length}; return 1;
}
static double delay_step(Delay *d,double input)
{
    /* 先读旧值再写新值，保证恰好延迟length拍。由delay_init初始化。 */
    double result=d->data[d->index]; d->data[d->index]=input;
    d->index=(d->index+1)%d->length; return result;
}
typedef struct { double tracking_iae,final,peak; } Result;
static Result run(int smith,int true_delay,double true_tau)
{
    Delay real,model;
    assert(delay_init(&real,true_delay) && delay_init(&model,100));
    PidState pid={0};
    const PidConfig cfg={1.5,.8,0,0,-3,3};
    const double dt=.01,a=exp(-dt/true_tau),am=exp(-dt);
    double y=0,free_model=0,delayed_model=0,area=0,peak=0;
    for(int k=0;k<4000;k++) {
        double feedback=smith ? free_model+y-delayed_model : y,u;
        assert(pid_update(&pid,&cfg,1,feedback,0,dt,&u));
        double real_input=delay_step(&real,u),model_input=delay_step(&model,u);
        double disturbance=k>=800 ? .2 : 0;
        y=a*y+(1-a)*(real_input+disturbance);
        free_model=am*free_model+(1-am)*u;
        delayed_model=am*delayed_model+(1-am)*model_input;
        if(k<800) { area+=fabs(1-y)*dt; peak=fmax(peak,y); }
        assert(isfinite(y) && fabs(y)<5);
        if(true_delay==100 && true_tau==1 && k<800)
            assert(fabs(y-delayed_model)<1e-12);
        /* 前100个区间真实对象仍没有收到任何输入。预测器不能抹去物理延迟。 */
        if(k<100) assert(y==0);
    }
    Result result={area,y,peak}; return result;
}
int main(void)
{
    Delay d;
    assert(delay_init(&d,3));
    assert(delay_step(&d,1)==0 && delay_step(&d,2)==0 && delay_step(&d,3)==0);
    assert(delay_step(&d,4)==1 && delay_step(&d,5)==2);
    assert(!delay_init(&d,0));
    Result plain=run(0,100,1),smith=run(1,100,1),mismatch=run(1,120,1.2);
    assert(smith.tracking_iae<plain.tracking_iae && smith.peak<plain.peak);
    assert(fabs(smith.final-1)<.001 && fabs(mismatch.final-1)<.005);
    printf("Smith: first 8s IAE PI=%.6f predictor=%.6f; peak=%.6f / %.6f\n",
           plain.tracking_iae,smith.tracking_iae,plain.peak,smith.peak);
    printf("Smith: final matched=%.6f, delay/tau mismatch=%.6f\n",smith.final,mismatch.final);
    return 0;
}

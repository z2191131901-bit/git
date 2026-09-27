/* 有限时长重复任务的模型逆ILC：每次试验后更新整段前馈。
 * 对象 y[k+1]=a*y[k]+b*u[k]+d[k]，u=feedforward+kp*(r-y)。
 * 每次试验重置物理初态，但保留学习前馈；不是连续周期重复控制器。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
enum { SAMPLES=400 };
static const double a=.95,b=.05,kp=2,learning=.6;
static double reference(int k)
{
    const double pi=3.14159265358979323846;
    return .5*(1-cos(2*pi*k/SAMPLES)); /* 起终点位置0，中间到1。 */
}
static double trial(const double feedforward[SAMPLES],double error[SAMPLES+1])
{
    double y=0,sum=0;
    for(int k=0;k<SAMPLES;k++) {
        error[k]=reference(k)-y;
        double input=feedforward[k]+kp*error[k];
        /* 与试验编号无关的重复扰动；换一次试验也必须保持相同时间对齐。 */
        double disturbance=.01*sin(4*3.14159265358979323846*k/SAMPLES);
        y=a*y+b*input+disturbance;
        assert(isfinite(y) && fabs(input)<10); /* 只检查范围，此处没有执行器限幅。 */
        double next_error=reference(k+1)-y; sum+=next_error*next_error;
    }
    error[SAMPLES]=reference(SAMPLES)-y;
    return sqrt(sum/SAMPLES);
}
static int learn(double feedforward[SAMPLES],const double error[SAMPLES+1],
                 double rate,double model_b)
{
    if(!feedforward || !error || !isfinite(rate) || rate<=0 || rate>=2 ||
       !isfinite(model_b) || model_b<=0) return 0;
    double next[SAMPLES];
    const double closed_a=a-b*kp;
    for(int k=0;k<SAMPLES;k++) {
        /* 相对阶为1：u[k]作用到y[k+1]，必须使用下一拍误差。
         * 已完成整次试验，读取error[k+1]不违反在线因果性。 */
        next[k]=feedforward[k]+rate*(error[k+1]-closed_a*error[k])/model_b;
        if(!isfinite(next[k])) return 0;
    }
    for(int k=0;k<SAMPLES;k++) feedforward[k]=next[k];
    return 1;
}
int main(void)
{
    double ff[SAMPLES]={0},error[SAMPLES+1],first,last,previous;
    first=previous=trial(ff,error); last=first;
    for(int iteration=1;iteration<=10;iteration++) {
        assert(learn(ff,error,learning,b));
        last=trial(ff,error);
        /* 精确模型、初态一致、重复扰动下，整条误差乘(1-learning)=.4。 */
        assert(fabs(last/previous-.4)<1e-8);
        previous=last;
    }
    assert(last<first*.001);
    double saved=ff[0];
    assert(!learn(ff,error,.6,0) && ff[0]==saved);
    /* 新开任务且清掉学习前馈，应回到最初表现。 */
    double empty[SAMPLES]={0},baseline=trial(empty,error);
    assert(fabs(baseline-first)<1e-12);
    printf("ILC: RMS initial=%.9f after 10 updates=%.9f; theoretical ratio=%.9f\n",
           first,last,pow(.4,10));
    return 0;
}

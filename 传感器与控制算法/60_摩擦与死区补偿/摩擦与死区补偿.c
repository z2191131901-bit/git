/* 简化一维电机：执行器输入死区 + 粘性/库仑摩擦 + 静摩擦。
 * 两组使用相同PI；补偿组加入参考速度前馈和执行器死区逆映射。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "../公共代码/PID核心.h"
static double deadzone(double command)
{
    return fabs(command)<=.2 ? 0 : copysign(fabs(command)-.2,command);
}
static double inverse_deadzone(double torque)
{
    return torque==0 ? 0 : torque+copysign(.2,torque);
}
static double plant(double velocity,double command,double dt)
{
    double force=deadzone(command),direction;
    /* 速度为0时，驱动力不足静摩擦阈值就保持不动。 */
    if(fabs(velocity)<1e-12) {
        if(fabs(force)<=.3) return 0;
        direction=copysign(1,force);
    } else direction=copysign(1,velocity);
    double next=velocity+dt*(force-velocity-.25*direction);
    /* 低驱动力下跨过0时捕获到静止，避免离散摩擦把速度反复推过零。 */
    if(velocity*next<0 && fabs(force)<=.3) next=0;
    return next;
}
static double run(int compensate,double *final)
{
    PidState s={0}; const PidConfig cfg={1,.8,0,0,-1,1};
    double v=0,area=0;
    const double reference=.1,dt=.001;
    for(int k=0;k<12000;k++) {
        /* tanh平滑摩擦方向；r=0时前馈为0，不凭噪声给任意方向的力。 */
        double ff=compensate ? reference+.25*tanh(reference/.02) : 0;
        double requested;
        assert(pid_update(&s,&cfg,reference,v,ff,dt,&requested));
        double command=compensate ? inverse_deadzone(requested) : requested;
        assert(fabs(command)<=1.2+1e-12);
        v=plant(v,command,dt);
        if(k<3000) area+=fabs(reference-v)*dt;
    }
    *final=v; return area;
}
int main(void)
{
    assert(deadzone(.2)==0 && deadzone(-.2)==0);
    assert(fabs(deadzone(.55)-.35)<1e-12);
    for(int i=-10;i<=10;i++) {
        double t=i*.05;
        assert(fabs(deadzone(inverse_deadzone(t))-t)<1e-12);
    }
    assert(plant(0,.49,.001)==0); /* 有效驱动力.29不足静摩擦.3。 */
    assert(plant(0,.6,.001)>0 && plant(0,-.6,.001)<0);
    assert(inverse_deadzone(0)==0);
    double end1,end2,area1=run(0,&end1),area2=run(1,&end2);
    assert(area2<area1*.5 && fabs(end2-.1)<.002);
    printf("friction compensation: first 3s IAE plain=%.6f compensated=%.6f; final=%.6f / %.6f\n",
           area1,area2,end1,end2);
    return 0;
}

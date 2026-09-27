/* 一阶直接MRAC：实际 y'=-a*y+b*u，参考 ym'=-am*ym+bm*r。
 * u=theta_y*y+theta_r*r；已知输入增益b为正。
 * 连续自适应律 theta'=-gamma*e*[y,r]，这里用小步长Euler离散。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
typedef struct { double theta_y,theta_r; } Adaptive;
static int adapt(Adaptive *s,double error,double y,double r,double gamma,double dt)
{
    if(!s || !isfinite(error) || !isfinite(y) || !isfinite(r) ||
       !isfinite(gamma) || gamma<0 || !isfinite(dt) || dt<=0) return 0;
    Adaptive next={s->theta_y-gamma*error*y*dt,s->theta_r-gamma*error*r*dt};
    if(!isfinite(next.theta_y) || !isfinite(next.theta_r)) return 0;
    *s=next; return 1;
}
typedef struct { double rms; Adaptive parameter; } Result;
static Result run(double gamma)
{
    const double dt=.001,plant_a=1,plant_b=.7,model_a=2,model_b=2;
    const double plant_decay=exp(-plant_a*dt),model_decay=exp(-model_a*dt);
    Adaptive s={0,1}; double y=0,ym=0,sum=0;
    for(int k=0;k<120000;k++) {
        double t=k*dt,r=.8*sin(.7*t)+.4*sin(1.7*t);
        double error=y-ym;
        double u=s.theta_y*y+s.theta_r*r;
        assert(isfinite(u) && fabs(u)<20);
        /* 先用旧参数发出控制，再依据同一时刻误差更新下一拍参数。 */
        assert(adapt(&s,error,y,r,gamma,dt));
        y=plant_decay*y+(1-plant_decay)*plant_b/plant_a*u;
        ym=model_decay*ym+(1-model_decay)*model_b/model_a*r;
        if(k>=100000) sum+=(y-ym)*(y-ym);
    }
    Result result={sqrt(sum/20000),s}; return result;
}
int main(void)
{
    Adaptive s={0,1};
    assert(adapt(&s,.2,.5,1,3,.01));
    assert(fabs(s.theta_y+.003)<1e-12 && fabs(s.theta_r-.994)<1e-12);
    assert(!adapt(&s,0,0,0,3,0) && fabs(s.theta_y+.003)<1e-12);
    Result fixed=run(0),adaptive=run(5);
    assert(adaptive.rms<.01 && adaptive.rms<fixed.rms*.1);
    /* 理想匹配参数theta_y=(a-am)/b，theta_r=bm/b；这里不宣称一般参数必收敛。 */
    printf("MRAC: last 20s RMS fixed=%.9f adaptive=%.9f\n",fixed.rms,adaptive.rms);
    printf("MRAC: theta_y=%.6f theta_r=%.6f; ideal=%.6f / %.6f\n",
           adaptive.parameter.theta_y,adaptive.parameter.theta_r,-1/.7,2/.7);
    return 0;
}

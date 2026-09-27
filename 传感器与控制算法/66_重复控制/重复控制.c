/* 连续运行的周期记忆控制：对象不在周期边界复位。
 * 已知闭环一阶模型的逆补偿，N拍记忆；与第65篇离线ILC更新不同。 */
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
enum { PERIOD=100 };
typedef struct {
    double command_memory[PERIOD];
    double error_history[PERIOD+1];
    size_t tick;
} Repetitive;
static int control(Repetitive *s,double error,double retention,int enabled,double *output)
{
    const double a=.9,b=.1,kp=1,learning=.4;
    if(!s || !output || !isfinite(error) || !isfinite(retention) ||
       retention<=0 || retention>1 || s->tick==SIZE_MAX) return 0;
    size_t k=s->tick,phase=k%PERIOD;
    double memory=0;
    if(enabled && k>=PERIOD) {
        /* 只使用过去的误差。N+1长度避免周期最后一拍误读已覆盖的数据。 */
        double old=s->error_history[(k-PERIOD)%(PERIOD+1)];
        double later=s->error_history[(k-PERIOD+1)%(PERIOD+1)];
        memory=retention*s->command_memory[phase]+learning/b*(later-(a-b*kp)*old);
    }
    double u=memory+kp*error;
    if(!isfinite(memory) || !isfinite(u)) return 0;
    s->command_memory[phase]=memory;
    s->error_history[k%(PERIOD+1)]=error;
    s->tick++;
    *output=u; return 1;
}
static double run(int enabled,double retention)
{
    Repetitive s={0};
    double y=0,sum=0;
    for(int k=0;k<50*PERIOD;k++) {
        double angle=2*3.14159265358979323846*k/PERIOD;
        double reference=sin(angle),error=reference-y,u;
        assert(control(&s,error,retention,enabled,&u));
        assert(fabs(u)<10);
        if(k>=49*PERIOD) sum+=error*error;
        /* 一拍离散模型，扰动与参考具有已知的整数周期。 */
        y=.9*y+.1*u+.02*cos(2*angle);
    }
    return sqrt(sum/PERIOD);
}
int main(void)
{
    Repetitive s={0}; double u=0;
    /* 编码历史为递增整数，独立检查整周期两端的索引，无需依赖闭环结果。 */
    for(int k=0;k<2*PERIOD;k++) {
        assert(control(&s,(double)k,1,1,&u));
        if(k>=PERIOD) {
            double old=(double)(k-PERIOD),later=old+1;
            assert(fabs(u-(k+4*(later-.8*old)))<1e-10);
        }
    }
    size_t saved=s.tick;
    assert(!control(&s,NAN,1,1,&u) && s.tick==saved);
    double baseline=run(0,1),ideal=run(1,1),leaky=run(1,.98);
    assert(baseline>.05 && ideal<1e-6 && leaky<baseline && leaky>ideal);
    printf("repetitive: final-cycle RMS feedback=%.9f memory=%.9f retention0.98=%.9f\n",
           baseline,ideal,leaky);
    return 0;
}

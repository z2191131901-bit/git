/* Welford在线均值/方差；RMS由 E[x²]=Var(x)+E[x]² 得到。
 * 不保存历史数组，每来一个样本只更新计数、均值、离差平方和。 */
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
typedef struct { uint32_t n; double mean,m2; } Stats;
static int add(Stats *s,double x)
{
    if(!s || !isfinite(x) || s->n==UINT32_MAX) return 0;
    Stats next=*s;
    next.n++;
    double delta=x-s->mean;
    next.mean+=delta/next.n;
    /* 第二个差值必须使用新均值。 */
    next.m2+=delta*(x-next.mean);
    if(!isfinite(next.mean) || !isfinite(next.m2) || next.m2<0) return 0;
    *s=next;
    return 1;
}
static int result(const Stats *s,double *variance,double *sample_variance,double *rms)
{
    if(!s || s->n<2 || !variance || !sample_variance || !rms) return 0;
    *variance=s->m2/s->n; /* 把已经采集的这些点当作完整总体。 */
    *sample_variance=s->m2/(s->n-1); /* 独立同分布抽样下的无偏估计。 */
    *rms=hypot(s->mean,sqrt(*variance));
    return 1;
}
int main(void)
{
    Stats s={0,0,0};
    double v=0,sv=0,rms=0;
    assert(!result(&s,&v,&sv,&rms));
    for(int i=1;i<=4;i++) assert(add(&s,i));
    assert(result(&s,&v,&sv,&rms));
    assert(s.mean==2.5 && v==1.25 && fabs(sv-5.0/3)<1e-12);
    assert(fabs(rms-sqrt(7.5))<1e-12);
    assert(!add(&s,NAN) && s.n==4);
    Stats offset={0,0,0};
    for(int i=1;i<=4;i++) assert(add(&offset,1e9+i));
    assert(result(&offset,&v,&sv,&rms) && fabs(v-1.25)<1e-12);
    printf("online stats: PASS, mean=2.5 variance=1.25 RMS=%.6f\n",sqrt(7.5));
    return 0;
}

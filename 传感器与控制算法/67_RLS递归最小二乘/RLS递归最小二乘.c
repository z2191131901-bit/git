/* 识别离散一阶对象 y[k+1]=a*y[k]+b*u[k]。
 * theta=[a,b]，回归向量phi=[y[k],u[k]]；lambda<1逐渐忘记旧数据。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
typedef struct { double theta[2],p[2][2]; } Rls;
static Rls initialize(double scale)
{
    Rls s={{0,0},{{scale,0},{0,scale}}}; return s;
}
static int update(Rls *s,double y,double u,double next_y,double forgetting)
{
    if(!s || !isfinite(y) || !isfinite(u) || !isfinite(next_y) ||
       !isfinite(forgetting) || forgetting<=0 || forgetting>1) return 0;
    double phi[2]={y,u},p_phi[2]={0};
    for(int i=0;i<2;i++)
        for(int j=0;j<2;j++) p_phi[i]+=s->p[i][j]*phi[j];
    double denominator=forgetting+phi[0]*p_phi[0]+phi[1]*p_phi[1];
    if(!isfinite(denominator) || denominator<=0) return 0;
    double residual=next_y-s->theta[0]*y-s->theta[1]*u;
    Rls next=*s;
    for(int i=0;i<2;i++) {
        next.theta[i]+=p_phi[i]/denominator*residual;
        if(!isfinite(next.theta[i])) return 0;
        for(int j=0;j<2;j++) {
            next.p[i][j]=(s->p[i][j]-p_phi[i]*p_phi[j]/denominator)/forgetting;
            if(!isfinite(next.p[i][j])) return 0;
        }
    }
    next.p[0][1]=next.p[1][0]=(next.p[0][1]+next.p[1][0])/2;
    /* 2x2对称矩阵正定检查；不把所有对角线为正当成充分条件。 */
    double determinant=next.p[0][0]*next.p[1][1]-next.p[0][1]*next.p[1][0];
    if(next.p[0][0]<=0 || next.p[1][1]<=0 || !isfinite(determinant) || determinant<=0)
        return 0;
    *s=next; return 1;
}
int main(void)
{
    Rls hand=initialize(1);
    assert(update(&hand,1,0,2,1));
    assert(hand.theta[0]==1 && hand.theta[1]==0 && hand.p[0][0]==.5);
    assert(!update(&hand,1,0,2,0) && hand.theta[0]==1);
    Rls fixed=initialize(100),tracking=initialize(100);
    double y=0;
    for(int k=0;k<2000;k++) {
        double a=k<1000 ? .92 : .8,b=k<1000 ? .08 : .15;
        double u=sin(.17*k)+.5*sin(.043*k);
        double next_y=a*y+b*u;
        assert(update(&fixed,y,u,next_y,1));
        assert(update(&tracking,y,u,next_y,.98));
        if(k==999) {
            assert(fabs(tracking.theta[0]-.92)<1e-5);
            assert(fabs(tracking.theta[1]-.08)<1e-5);
        }
        y=next_y;
    }
    double fixed_error=hypot(fixed.theta[0]-.8,fixed.theta[1]-.15);
    double tracking_error=hypot(tracking.theta[0]-.8,tracking.theta[1]-.15);
    assert(tracking_error<1e-5 && tracking_error<fixed_error*.01);
    /* 没有激励时，不会凭空识别出参数；遗忘反而使P增大。 */
    Rls silent=initialize(1);
    for(int k=0;k<20;k++) assert(update(&silent,0,0,0,.98));
    assert(silent.theta[0]==0 && silent.theta[1]==0 && silent.p[0][0]>1);
    printf("RLS: final no forgetting a=%.6f b=%.6f; forgetting a=%.6f b=%.6f\n",
           fixed.theta[0],fixed.theta[1],tracking.theta[0],tracking.theta[1]);
    printf("RLS: parameter error=%.9f / %.9f; unexcited P00=%.6f\n",
           fixed_error,tracking_error,silent.p[0][0]);
    return 0;
}

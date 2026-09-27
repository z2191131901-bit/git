/* 双积分对象增加误差积分：z_next=z+dt*(r-p)。
 * 状态为[p-r,v,z]，设计3维离散LQR即得到LQI；无输入限幅。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
enum { MAX_N=3 };
typedef struct { double p[3][3],k[3]; int iterations; } Design;
static int design(int n,double dt,double r,Design *out)
{
    if(!out || (n!=2 && n!=3) || !isfinite(dt) || dt<=0 ||
       !isfinite(r) || r<=0) return 0;
    const double a[3][3]={{1,dt,0},{0,1,0},{-dt,0,1}};
    const double b[3]={dt*dt/2,dt,0},q[3]={10,1,5};
    Design d={0};
    for(int i=0;i<n;i++) d.p[i][i]=q[i];
    int converged=0;
    for(int iteration=0;iteration<20000;iteration++) {
        double pb[3]={0},bpa[3]={0},next[3][3]={{0}},denominator=r;
        for(int i=0;i<n;i++)
            for(int j=0;j<n;j++) pb[i]+=d.p[i][j]*b[j];
        for(int i=0;i<n;i++) denominator+=b[i]*pb[i];
        if(!isfinite(denominator) || denominator<=0) return 0;
        for(int j=0;j<n;j++)
            for(int i=0;i<n;i++) bpa[j]+=pb[i]*a[i][j];
        double change=0,scale=1;
        for(int i=0;i<n;i++)
            for(int j=0;j<n;j++) {
                for(int row=0;row<n;row++)
                    for(int col=0;col<n;col++)
                        next[i][j]+=a[row][i]*d.p[row][col]*a[col][j];
                next[i][j]-=bpa[i]*bpa[j]/denominator;
                if(i==j) next[i][j]+=q[i];
                if(!isfinite(next[i][j])) return 0;
                change=fmax(change,fabs(next[i][j]-d.p[i][j]));
                scale=fmax(scale,fabs(next[i][j]));
            }
        memcpy(d.p,next,sizeof next); d.iterations=iteration+1;
        if(change<1e-12*scale) { converged=1; break; }
    }
    if(!converged) return 0;
    /* 用最终P重新算K，不保留上一轮尚未收敛的增益。 */
    double pb[3]={0},denominator=r;
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++) pb[i]+=d.p[i][j]*b[j];
    for(int i=0;i<n;i++) denominator+=b[i]*pb[i];
    for(int j=0;j<n;j++) {
        for(int i=0;i<n;i++) d.k[j]+=pb[i]*a[i][j];
        d.k[j]/=denominator;
        if(!isfinite(d.k[j])) return 0;
    }
    *out=d; return 1;
}
static double value(const Design *d,const double x[3])
{
    double sum=0;
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) sum+=x[i]*d->p[i][j]*x[j];
    return sum;
}
static double run(const Design *d,int integral,double *last_u)
{
    const double dt=.05,reference=1;
    double position=0,velocity=0,z=0,u=0;
    for(int k=0;k<1600;k++) {
        u=-d->k[0]*(position-reference)-d->k[1]*velocity-d->k[2]*z;
        double acceleration=u+(k>=100 ? .3 : 0);
        /* 积分使用旧位置，与设计的增广A保持一致。 */
        if(integral) z+=dt*(reference-position);
        position+=velocity*dt+acceleration*dt*dt/2;
        velocity+=acceleration*dt;
        assert(isfinite(position) && fabs(position)<10);
    }
    *last_u=u; return position-reference;
}
int main(void)
{
    Design lqr,lqi;
    assert(design(2,.05,.1,&lqr) && design(3,.05,.1,&lqi));
    assert(lqi.k[2]<0);
    /* Bellman恒等式核对P、K和积分方向；要求残差接近0。 */
    double x[3]={.4,-.2,.3},u=0;
    for(int i=0;i<3;i++) u-=lqi.k[i]*x[i];
    double next[3]={x[0]+.05*x[1]+.05*.05*u/2,x[1]+.05*u,x[2]-.05*x[0]};
    double cost=10*x[0]*x[0]+x[1]*x[1]+5*x[2]*x[2]+.1*u*u;
    assert(fabs(value(&lqi,x)-value(&lqi,next)-cost)<1e-7);
    double u1,u2,e1=run(&lqr,0,&u1),e2=run(&lqi,1,&u2);
    assert(fabs(e1-.3/lqr.k[0])<1e-8);
    assert(fabs(e2)<1e-8 && fabs(u2+.3)<1e-8);
    double saved=lqi.k[0];
    assert(!design(3,.05,0,&lqi) && lqi.k[0]==saved);
    printf("LQI: K=[%.6f, %.6f, %.6f], iterations=%d\n",
           lqi.k[0],lqi.k[1],lqi.k[2],lqi.iterations);
    printf("LQI: disturbance position error LQR=%.6f LQI=%.9f; steady u=%.6f\n",e1,e2,u2);
    return 0;
}

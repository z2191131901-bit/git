/*
 * 离散 LQR：二维位置-速度状态、标量加速度输入的双积分模型。
 * x_next=A*x+B*u，A=[[1,dt],[0,1]]，B=[dt²/2,dt]^T。
 * 求无限时域离散 Riccati 方程，获得固定反馈 u=-K*x。
 * 此例没有输入约束，不能把限幅后的控制仍称为无约束最优解。
 */
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

typedef struct { double p[2][2], k[2]; int iterations; } LqrDesign;

/* 计算一次 Riccati 递推，同时返回当前 P 对应的反馈增益。 */
static int riccati_step(double dt, double q_position, double q_velocity, double r,
                        const double p[2][2], double next[2][2], double k[2])
{
    const double a[2][2]={{1,dt},{0,1}}, b[2]={dt*dt/2,dt};
    double pb[2]={p[0][0]*b[0]+p[0][1]*b[1],
                  p[1][0]*b[0]+p[1][1]*b[1]};
    const double denominator=r+b[0]*pb[0]+b[1]*pb[1];
    if (!isfinite(denominator) || denominator<=0) return 0;
    k[0]=pb[0]/denominator;
    k[1]=(pb[0]*dt+pb[1])/denominator;
    if (!isfinite(k[0]) || !isfinite(k[1])) return 0;

    for (int row=0;row<2;++row)
    {
        for (int col=0;col<2;++col)
        {
            /* A^T*P*A 的第 row,col 元素，四个索引对应矩阵乘法。 */
            double value=0;
            for (int i=0;i<2;++i)
                for (int j=0;j<2;++j)
                    value+=a[i][row]*p[i][j]*a[j][col];
            value-=denominator*k[row]*k[col];
            if (row==col) value+=row==0?q_position:q_velocity;
            if (!isfinite(value)) return 0;
            next[row][col]=value;
        }
    }
    /* 理论上 P 对称，消除浮点运算造成的极小不对称。 */
    next[0][1]=next[1][0]=(next[0][1]+next[1][0])/2;
    return 1;
}
static int lqr_design(double dt, double q_position, double q_velocity,
                      double r, LqrDesign *output)
{
    if (output==NULL || !isfinite(dt) || dt<=0 ||
        !isfinite(q_position) || q_position<=0 ||
        !isfinite(q_velocity) || q_velocity<0 || !isfinite(r) || r<=0) return 0;
    LqrDesign result={{{q_position,0},{0,q_velocity}},{0,0},0};
    int converged=0;
    for (int iteration=0;iteration<20000;++iteration)
    {
        double next[2][2],gain[2],change=0,scale=1;
        /* C99 的二维数组 const 转换使用显式指针类型，避免旧编译器诊断。 */
        int ok=riccati_step(dt,q_position,q_velocity,r,
                            (const double (*)[2])result.p,next,gain);
        if (!ok) return 0;
        for (int i=0;i<2;++i)
            for (int j=0;j<2;++j)
            {
                change=fmax(change,fabs(next[i][j]-result.p[i][j]));
                scale=fmax(scale,fabs(next[i][j]));
                result.p[i][j]=next[i][j];
            }
        result.iterations=iteration+1;
        if (change<1e-12*scale) { converged=1; break; }
    }
    if (!converged) return 0; /* 有迭代上限，不把未收敛结果当有效控制器。 */
    double unused[2][2];
    if (!riccati_step(dt,q_position,q_velocity,r,
                      (const double (*)[2])result.p,unused,result.k)) return 0;
    *output=result;
    return 1;
}
static double value_function(const LqrDesign *design, double position, double velocity)
{
    return design->p[0][0]*position*position
        +2*design->p[0][1]*position*velocity+design->p[1][1]*velocity*velocity;
}
int main(void)
{
    const double dt=0.05,qp=10,qv=1,r=0.1;
    LqrDesign design;
    int ok=lqr_design(dt,qp,qv,r,&design);
    assert(ok);
    /* 检查闭环 A-BK 的两个特征值均在单位圆内。 */
    const double a00=1-dt*dt/2*design.k[0],a01=dt-dt*dt/2*design.k[1];
    const double a10=-dt*design.k[0],a11=1-dt*design.k[1];
    const double trace=a00+a11,determinant=a00*a11-a01*a10;
    const double discriminant=trace*trace-4*determinant;
    const double radius=discriminant<0?sqrt(determinant):
        fmax(fabs((trace+sqrt(discriminant))/2),fabs((trace-sqrt(discriminant))/2));
    assert(radius<1);

    double position=1,velocity=-0.2;
    const double u=-design.k[0]*position-design.k[1]*velocity;
    const double next_position=position+dt*velocity+dt*dt/2*u;
    const double next_velocity=velocity+dt*u;
    const double bellman_residual=value_function(&design,position,velocity)
        -(qp*position*position+qv*velocity*velocity+r*u*u
          +value_function(&design,next_position,next_velocity));
    assert(fabs(bellman_residual)<1e-7);

    for (int i=0;i<400;++i)
    {
        const double control=-design.k[0]*position-design.k[1]*velocity;
        position+=dt*velocity+dt*dt/2*control;
        velocity+=dt*control;
    }
    assert(fabs(position)<1e-6 && fabs(velocity)<1e-6);
    double saved=design.k[0];
    ok=lqr_design(dt,qp,qv,0,&design);
    assert(!ok && design.k[0]==saved);
    printf("LQR K=[%.6f, %.6f], spectral radius=%.6f, iterations=%d\n",
           design.k[0],design.k[1],radius,design.iterations);
    puts("lqr: PASS");
    return 0;
}

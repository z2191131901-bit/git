/*
 * 两步连续输入约束 MPC：
 * x_next=0.9*x+0.2*u，0<=u<=1，0<=预测状态<=1.2。
 * 优化变量为 u0/u1；精确枚举二维凸二次规划的内部、边界和顶点候选。
 * 不是控制量网格搜索，也不是可直接替代通用求解器的工业 MPC 库。
 */
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

static const double MODEL_A=0.9,MODEL_B=0.2,CONTROL_WEIGHT=0.02;
static const double U_MIN=0,U_MAX=1,X_MIN=0,X_MAX=1.2;
typedef struct { double c0,c1,bound; } Constraint; /* c0*u0+c1*u1<=bound。 */
typedef struct
{
    double u0,u1;        /* 两步计划，只执行u0，下个采样重新求解。 */
    double x1,x2,cost;
} MpcPlan;

static void make_constraints(double state, Constraint constraints[8])
{
    constraints[0]=(Constraint){1,0,U_MAX};
    constraints[1]=(Constraint){-1,0,-U_MIN};
    constraints[2]=(Constraint){0,1,U_MAX};
    constraints[3]=(Constraint){0,-1,-U_MIN};
    /* x1=A*x+B*u0，x2=A²*x+A*B*u0+B*u1。 */
    constraints[4]=(Constraint){MODEL_B,0,X_MAX-MODEL_A*state};
    constraints[5]=(Constraint){-MODEL_B,0,MODEL_A*state-X_MIN};
    constraints[6]=(Constraint){MODEL_A*MODEL_B,MODEL_B,X_MAX-MODEL_A*MODEL_A*state};
    constraints[7]=(Constraint){-MODEL_A*MODEL_B,-MODEL_B,MODEL_A*MODEL_A*state-X_MIN};
}
static double plan_cost(double state, double reference, double u0, double u1)
{
    const double equilibrium=(1-MODEL_A)*reference/MODEL_B;
    const double x1=MODEL_A*state+MODEL_B*u0;
    const double x2=MODEL_A*x1+MODEL_B*u1;
    /*
     * 惩罚对参考平衡输入的偏离，而不是一直把输入往0拉。
     * 对本模型，达到非零参考本来就需要非零维持输入。
     */
    return (x1-reference)*(x1-reference)+(x2-reference)*(x2-reference)
        +CONTROL_WEIGHT*((u0-equilibrium)*(u0-equilibrium)
                         +(u1-equilibrium)*(u1-equilibrium));
}
static void consider(double state, double reference, const Constraint constraints[8],
                      double u0, double u1, MpcPlan *best)
{
    if (!isfinite(u0) || !isfinite(u1)) return;
    for (int i=0;i<8;++i)
        if (constraints[i].c0*u0+constraints[i].c1*u1>constraints[i].bound+1e-9)
            return;
    const double cost=plan_cost(state,reference,u0,u1);
    if (!isfinite(cost) || cost>=best->cost) return;
    *best=(MpcPlan){u0,u1,MODEL_A*state+MODEL_B*u0,
                   MODEL_A*MODEL_A*state+MODEL_A*MODEL_B*u0+MODEL_B*u1,cost};
}
static int mpc_solve(double state, double reference, MpcPlan *output)
{
    if (output==NULL || !isfinite(state) || !isfinite(reference) ||
        reference<X_MIN || reference>X_MAX) return 0;
    const double equilibrium=(1-MODEL_A)*reference/MODEL_B;
    if (equilibrium<U_MIN || equilibrium>U_MAX) return 0;
    Constraint constraints[8];
    make_constraints(state,constraints);
    /*
     * J = 0.5*U^T*H*U + g^T*U + 与U无关的常数。
     * 预测矩阵为 [[B,0],[A*B,B]]；H=2*(G^T*G+R*I)。
     */
    const double b=MODEL_B,ab=MODEL_A*MODEL_B;
    const double h00=2*(b*b+ab*ab+CONTROL_WEIGHT),h01=2*ab*b;
    const double h11=2*(b*b+CONTROL_WEIGHT);
    const double d1=MODEL_A*state-reference,d2=MODEL_A*MODEL_A*state-reference;
    const double g0=2*(b*d1+ab*d2-CONTROL_WEIGHT*equilibrium);
    const double g1=2*(b*d2-CONTROL_WEIGHT*equilibrium);
    const double determinant=h00*h11-h01*h01;
    const double inverse00=h11/determinant,inverse01=-h01/determinant;
    const double inverse11=h00/determinant;
    const double free0=-(inverse00*g0+inverse01*g1);
    const double free1=-(inverse01*g0+inverse11*g1);
    MpcPlan best={0,0,0,0,INFINITY};

    /* ① 内部候选：无约束极小值 -H^-1*g。 */
    consider(state,reference,constraints,free0,free1,&best);
    for (int i=0;i<8;++i)
    {
        Constraint c=constraints[i];
        /*
         * ② 一条边界上的极小值：
         * U=U_free-H^-1*c*lambda
         * lambda=(c^T*U_free-bound)/(c^T*H^-1*c)。
         */
        const double v0=inverse00*c.c0+inverse01*c.c1;
        const double v1=inverse01*c.c0+inverse11*c.c1;
        const double denominator=c.c0*v0+c.c1*v1;
        const double lambda=(c.c0*free0+c.c1*free1-c.bound)/denominator;
        consider(state,reference,constraints,free0-v0*lambda,free1-v1*lambda,&best);

        /* ③ 两条非平行边界交点，是可行多边形的潜在顶点。 */
        for (int j=i+1;j<8;++j)
        {
            Constraint d=constraints[j];
            const double det=c.c0*d.c1-c.c1*d.c0;
            if (fabs(det)<1e-12) continue;
            const double u0=(c.bound*d.c1-c.c1*d.bound)/det;
            const double u1=(c.c0*d.bound-c.bound*d.c0)/det;
            consider(state,reference,constraints,u0,u1,&best);
        }
    }
    if (!isfinite(best.cost)) return 0; /* 无可行方案，不修改旧输出。 */
    *output=best;
    return 1;
}
/* 独立网格仅用于检查：解析求解的代价不应比任何可行网格方案更差。 */
static void check_against_grid(double state, double reference, const MpcPlan *plan)
{
    for (int i=0;i<=100;++i)
        for (int j=0;j<=100;++j)
        {
            const double u0=i/100.0,u1=j/100.0;
            const double x1=0.9*state+0.2*u0,x2=0.9*x1+0.2*u1;
            if (x1<0 || x1>1.2 || x2<0 || x2>1.2) continue;
            assert(plan->cost<=plan_cost(state,reference,u0,u1)+1e-9);
        }
}
int main(void)
{
    MpcPlan plan;
    int ok=mpc_solve(0.8,0.8,&plan);
    assert(ok && fabs(plan.u0-0.4)<1e-12 && fabs(plan.u1-0.4)<1e-12);
    assert(fabs(plan.cost)<1e-12);
    check_against_grid(0.8,0.8,&plan);

    ok=mpc_solve(0,0.8,&plan);
    assert(ok && fabs(plan.u0-1)<1e-12 && fabs(plan.u1-1)<1e-12);
    assert(fabs(plan.cost-0.5508)<1e-12);
    check_against_grid(0,0.8,&plan);

    /* 当前状态1.3虽在范围外，但一步可恢复；x1上界约束会激活。 */
    ok=mpc_solve(1.3,1.2,&plan);
    assert(ok && fabs(plan.u0-0.15)<1e-10 && fabs(plan.x1-1.2)<1e-10);
    assert(fabs(plan.cost-0.00405)<1e-12);
    check_against_grid(1.3,1.2,&plan);
    printf("state-constrained plan: u0=%.6f u1=%.6f x1=%.6f\n",plan.u0,plan.u1,plan.x1);

    MpcPlan saved=plan;
    ok=mpc_solve(1.5,0.8,&plan);
    assert(!ok && plan.u0==saved.u0); /* 即使u=0，x1=1.35仍超上界，无解。 */

    double state=0;
    for (int sample=0;sample<100;++sample)
    {
        ok=mpc_solve(state,0.8,&plan);
        assert(ok && plan.u0>=-1e-9 && plan.u0<=1+1e-9);
        state=0.9*state+0.2*plan.u0; /* 只执行第一步，随后重新测量和优化。 */
        assert(state>=-1e-9 && state<=1.2+1e-9);
    }
    assert(fabs(state-0.8)<1e-8);
    printf("MPC closed-loop final state=%.8f reference=0.8\n",state);
    puts("mpc: PASS");
    return 0;
}

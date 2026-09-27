/*
 * 五状态平面 EKF：s=[x,y,yaw,v,bg]。
 * 预测输入：Z 轴角速度、已校准且去重力的前向加速度。
 * 可选观测：编码器前向速度、磁力计处理后的航向、外部局部 x/y。
 * 本例限定水平无侧滑车辆；不是完整三维 INS/ESKF。
 * 每个数据到达时调用对应更新，没数据就跳过，不拿 0 冒充观测。
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
enum { PX, PY, YAW, SPEED, BIAS, N };
typedef struct { double s[N], p[N][N]; } Filter;
static double wrap(double a) { return atan2(sin(a),cos(a)); }

static int valid(const Filter *f)
{
    if (!f) return 0;
    for (int i=0;i<N;i++) {
        if (!isfinite(f->s[i]) || !isfinite(f->p[i][i]) || f->p[i][i]<0) return 0;
        for (int j=0;j<N;j++)
            if (!isfinite(f->p[i][j]) || fabs(f->p[i][j]-f->p[j][i])>1e-9)
                return 0;
    }
    /* 调用方还必须提供半正定初始 P；有限性/对称性不能证明半正定。 */
    return 1;
}
static Filter initialize(void)
{
    Filter f={0};
    const double variance[N]={1,1,.25,.25,.01};
    for (int i=0;i<N;i++) f.p[i][i]=variance[i];
    return f;
}
static int predict(Filter *f, double gyro_z, double forward_accel, double dt)
{
    if (!valid(f) || !isfinite(gyro_z) || !isfinite(forward_accel) ||
        !isfinite(dt) || dt<=0 || dt>0.1) return 0;
    Filter next=*f;
    double transition[N][N]={{0}}, product[N][N]={{0}};
    const double yaw=f->s[YAW], v=f->s[SPEED];
    const double c=cos(yaw), s=sin(yaw);
    /* 为看清雅可比，使用一阶欧拉递推。减小 dt 可降低离散误差。 */
    next.s[PX]+=v*c*dt;
    next.s[PY]+=v*s*dt;
    next.s[YAW]=wrap(yaw+(gyro_z-f->s[BIAS])*dt);
    next.s[SPEED]+=forward_accel*dt;
    /* bg 采用随机游走模型：均值不变，不确定度通过 Q 增长。 */
    for (int i=0;i<N;i++) transition[i][i]=1;
    transition[PX][YAW]=-v*s*dt;
    transition[PX][SPEED]=c*dt;
    transition[PY][YAW]=v*c*dt;
    transition[PY][SPEED]=s*dt;
    transition[YAW][BIAS]=-dt;
    /* P_pred = F*P*F^T + Q；不能只更新 P 的对角线。 */
    for (int i=0;i<N;i++)
        for (int j=0;j<N;j++)
            for (int k=0;k<N;k++)
                product[i][j]+=transition[i][k]*f->p[k][j];
    for (int i=0;i<N;i++)
        for (int j=0;j<N;j++) {
            next.p[i][j]=0;
            for (int k=0;k<N;k++)
                next.p[i][j]+=product[i][k]*transition[j][k];
        }
    /* 对角 Q 是教学近似，数值为每秒增长率，单位见讲解。
     * 不把它声称为严格连续白噪声模型离散化得到的 Q。 */
    const double q_rate[N]={.0025,.0025,.0004,.04,.000001};
    for (int i=0;i<N;i++) next.p[i][i]+=q_rate[i]*dt;
    if (!valid(&next)) return 0;
    *f=next;
    return 1;
}
/*
 * 一次观测一个状态分量，H 是只有 index 项为 1 的行向量。
 * return 1=接受，0=异常值门控拒绝，-1=参数错误；拒绝/错误不改状态。
 * variance 是测量方差，不是标准差；gate2=9 表示标量 3 sigma 门限。
 */
static int observe(Filter *f, int index, double measurement,
                   double variance, double gate2)
{
    if (!valid(f) || index<0 || index>=N || !isfinite(measurement) ||
        !isfinite(variance) || variance<=0 || !isfinite(gate2) || gate2<=0)
        return -1;
    double residual=measurement-f->s[index];
    if (index==YAW) residual=wrap(residual); /* 179° 与 -179° 只差 2°。 */
    const double innovation_variance=f->p[index][index]+variance;
    if (!isfinite(residual) || !isfinite(innovation_variance) ||
        innovation_variance<=0) return -1;
    if (fabs(residual)/sqrt(innovation_variance)>sqrt(gate2)) return 0;
    Filter next=*f;
    double gain[N], a[N][N]={{0}}, ap[N][N]={{0}};
    for (int i=0;i<N;i++) {
        gain[i]=f->p[i][index]/innovation_variance;
        next.s[i]+=gain[i]*residual;
        a[i][i]=1;
        a[i][index]-=gain[i]; /* A = I-KH。 */
    }
    next.s[YAW]=wrap(next.s[YAW]);
    /* Joseph 形式：P_new=(I-KH)P(I-KH)^T + K R K^T。
     * 小矩阵多算几次乘法，换取更可靠的协方差数值表现。 */
    for (int i=0;i<N;i++)
        for (int j=0;j<N;j++)
            for (int k=0;k<N;k++) ap[i][j]+=a[i][k]*f->p[k][j];
    for (int i=0;i<N;i++)
        for (int j=0;j<N;j++) {
            next.p[i][j]=gain[i]*variance*gain[j];
            for (int k=0;k<N;k++) next.p[i][j]+=ap[i][k]*a[j][k];
        }
    /* 去掉浮点运算造成的极小非对称误差。 */
    for (int i=0;i<N;i++)
        for (int j=i+1;j<N;j++)
            next.p[i][j]=next.p[j][i]=.5*(next.p[i][j]+next.p[j][i]);
    if (!valid(&next)) return -1;
    *f=next;
    return 1;
}
/* 测试用 Cholesky：检查整个矩阵正定，不能只看对角线是否非负。 */
static void check_covariance(const Filter *f)
{
    double l[N][N]={{0}};
    assert(valid(f));
    for (int i=0;i<N;i++)
        for (int j=0;j<=i;j++) {
            double value=f->p[i][j];
            for (int k=0;k<j;k++) value-=l[i][k]*l[j][k];
            if (i==j) { assert(value>0); l[i][j]=sqrt(value); }
            else l[i][j]=value/l[j][j];
        }
}
int main(void)
{
    /* 手算：位置先验方差 1，测量方差 1，观测 2，K=.5，结果 x=1、P=.5。 */
    Filter unit=initialize();
    assert(observe(&unit,PX,2,1,9)==1);
    assert(fabs(unit.s[PX]-1)<1e-12 && fabs(unit.p[PX][PX]-.5)<1e-12);
    Filter saved=unit;
    assert(observe(&unit,PX,1000,1,9)==0);
    assert(memcmp(&unit,&saved,sizeof unit)==0);
    assert(observe(&unit,PX,0,0,9)==-1);
    assert(!predict(&unit,0,0,0));
    assert(memcmp(&unit,&saved,sizeof unit)==0);
    const double pi=3.14159265358979323846;
    unit.s[YAW]=179*pi/180;
    assert(observe(&unit,YAW,-179*pi/180,.01,9)==1);
    assert(fabs(wrap(unit.s[YAW]-pi))<2*pi/180);

    Filter fused=initialize(), dead=initialize();
    fused.s[SPEED]=dead.s[SPEED]=1;
    /* 独立真值用解析圆周公式生成，避免拿滤波器递推式充当真值。 */
    const double dt=.01, true_rate=.08, true_bias=.02;
    double truth_x=0, truth_y=0, truth_yaw=0;
    double outage_error=0, recovery_error=0;
    int rejected=0, positions=0;
    for (int k=1;k<=3000;k++) {
        const double t=k*dt;
        truth_yaw=wrap(true_rate*t);
        truth_x=sin(true_rate*t)/true_rate;
        truth_y=(1-cos(true_rate*t))/true_rate;
        /* 固定正弦扰动便于每次复现；不是统计蒙特卡洛实验。 */
        const double gyro=true_rate+true_bias+.001*sin(3*t);
        const double accel=.005*sin(2*t);
        assert(predict(&fused,gyro,accel,dt));
        assert(predict(&dead,gyro,accel,dt));
        /* 轮速 50 Hz：融合速度，别再把同一轮速积分的位置当独立观测。 */
        if (k%2==0) {
            const double wheel=1+.01*sin(4*t);
            assert(observe(&fused,SPEED,wheel,.02*.02,9)==1);
            assert(observe(&dead,SPEED,wheel,.02*.02,9)==1);
        }
        /* 磁航向 10 Hz；第 8 秒故意制造约 2 rad 的磁干扰。 */
        if (k%10==0) {
            const double magnetic=wrap(truth_yaw+.01*sin(t)+(k==800 ? 2 : 0));
            const int accepted=observe(&fused,YAW,magnetic,.03*.03,9);
            assert(accepted>=0);
            if (k==800) { assert(accepted==0); rejected++; }
        }
        /* 外部位置 5 Hz：模拟 GNSS/UWB/定位模块已输出同一坐标系下米制 x/y。
         * 10~18 秒完全断开；第 24 秒加 100 m 跳变，分别进行标量门控。 */
        if (k%20==0 && !(k>=1000 && k<=1800)) {
            const double jump=k==2400 ? 100 : 0;
            int ax=observe(&fused,PX,truth_x+.03*sin(2*t)+jump,.1*.1,9);
            int ay=observe(&fused,PY,truth_y+.03*cos(2*t)+jump,.1*.1,9);
            assert(ax>=0 && ay>=0);
            if (k==2400) { assert(ax==0 && ay==0); rejected+=2; }
            else { assert(ax==1 && ay==1); positions++; }
        }
        if (k==1800) outage_error=hypot(fused.s[PX]-truth_x,fused.s[PY]-truth_y);
        if (k==2200) recovery_error=hypot(fused.s[PX]-truth_x,fused.s[PY]-truth_y);
        check_covariance(&fused);
        check_covariance(&dead);
    }
    const double fused_error=hypot(fused.s[PX]-truth_x,fused.s[PY]-truth_y);
    const double dead_error=hypot(dead.s[PX]-truth_x,dead.s[PY]-truth_y);
    assert(rejected==3 && positions>0);
    assert(fused_error<.15 && fused_error<dead_error*.1);
    assert(fabs(fused.s[BIAS]-true_bias)<.005);
    assert(outage_error<.5 && recovery_error<.15);
    printf("navigation EKF: dead=%.6f m fused=%.6f m bias=%.6f rad/s\n",
           dead_error,fused_error,fused.s[BIAS]);
    printf("outage end=%.6f m recovery=%.6f m rejected scalar observations=%d\n",
           outage_error,recovery_error,rejected);
    return 0;
}

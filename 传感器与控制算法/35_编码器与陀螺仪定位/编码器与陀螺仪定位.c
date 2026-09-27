/*
 * 编码器距离 + 陀螺仪转角：平面航位推算。
 * 机体 X 向前、Y 向左；世界 X/Y 为初始化时的局部坐标。
 * ds 单位 m，gyro_z 单位 rad/s，dt 单位 s，yaw 单位 rad。
 * 编码器值是两个轮子在同一个采样区间内的有符号行程，倒车为负。
 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
typedef struct { double x, y, yaw; } Pose2;
static double wrap(double a) { return atan2(sin(a),cos(a)); }

/* 返回 0 表示输入不合法，原状态保持不变。bias 是预先标定的 Z 轴零偏。 */
static int odometry_step(Pose2 *state, double left_m, double right_m,
                         double gyro_z, double bias, double dt)
{
    if (!state || !isfinite(left_m) || !isfinite(right_m) ||
        !isfinite(gyro_z) || !isfinite(bias) || !isfinite(dt) || dt<=0 ||
        !isfinite(state->x) || !isfinite(state->y) || !isfinite(state->yaw))
        return 0;
    const double ds=0.5*(left_m+right_m);
    const double angle=(gyro_z-bias)*dt;
    const double half=0.5*angle;
    /* 恒定线速度/角速度时，弧长转换成弦长：
     * chord = ds * sin(angle/2)/(angle/2)。
     * 直接用 ds*cos(yaw+angle/2) 是中点近似，这里补上 sinc 因子。 */
    const double sinc=fabs(half)<1e-6 ? 1-half*half/6 : sin(half)/half;
    Pose2 next=*state;
    next.x+=ds*sinc*cos(state->yaw+half);
    next.y+=ds*sinc*sin(state->yaw+half);
    next.yaw=wrap(state->yaw+angle);
    if (!isfinite(next.x) || !isfinite(next.y) || !isfinite(next.yaw))
        return 0;
    *state=next;
    return 1;
}
int main(void)
{
    const double pi=3.14159265358979323846;
    Pose2 straight={0,0,0}, reverse={0,0,0}, circle={0,0,0};
    for (int i=0;i<100;i++)
        assert(odometry_step(&straight,0.01,0.01,0.02,0.02,0.01));
    assert(fabs(straight.x-1)<1e-12 && fabs(straight.y)<1e-12);
    assert(odometry_step(&reverse,-1,-1,0,0,1));
    assert(fabs(reverse.x+1)<1e-12);
    /* v=1 m/s，omega=1 rad/s，走 pi/2 秒，理论终点 (1,1)。 */
    for (int i=0;i<100;i++) {
        double dt=pi/200;
        assert(odometry_step(&circle,dt,dt,1.02,0.02,dt));
    }
    assert(fabs(circle.x-1)<1e-12 && fabs(circle.y-1)<1e-12);
    Pose2 spin={0,0,0};
    assert(odometry_step(&spin,-0.1,0.1,1,0,1));
    assert(spin.x==0 && spin.y==0 && fabs(spin.yaw-1)<1e-12);
    Pose2 saved=spin;
    assert(!odometry_step(&spin,NAN,0,0,0,0.01));
    assert(spin.x==saved.x && spin.yaw==saved.yaw);
    assert(!odometry_step(&spin,0,0,0,0,0));
    /* 直线走 10 秒，故意保留 0.02 rad/s 零偏，展示侧向漂移。 */
    Pose2 uncorrected={0,0,0}, corrected={0,0,0};
    for (int i=0;i<1000;i++) {
        assert(odometry_step(&uncorrected,.01,.01,.02,0,.01));
        assert(odometry_step(&corrected,.01,.01,.02,.02,.01));
    }
    assert(fabs(uncorrected.y)>.9 && fabs(corrected.y)<1e-12);
    printf("encoder+gyro: quarter circle=(%.6f, %.6f); bias drift y=%.6f, calibrated y=%.6f\n",
           circle.x,circle.y,uncorrected.y,corrected.y);
    return 0;
}

/*
 * 差速小车里程计：左右轮走过的距离 -> 平面位置和朝向。
 * 世界 X 为初始前方、Y 为左方，yaw 逆时针为正；距离 m、角度 rad。
 * 假设刚性平面、无侧滑，采样间隔内左右轮速度比近似不变。
 */
#include "../common/math3d.h"
#include <assert.h>
#include <stdio.h>

typedef struct { double x,y,yaw; } Pose2D;
static int odometry_update(Pose2D *pose, double left_distance,
                            double right_distance, double track_width)
{
    if (pose==NULL || !isfinite(pose->x) || !isfinite(pose->y) || !isfinite(pose->yaw) ||
        !isfinite(left_distance) || !isfinite(right_distance) ||
        !isfinite(track_width) || track_width<=0) return 0;
    const double distance=left_distance/2+right_distance/2;
    const double angle=(right_distance-left_distance)/track_width;
    if (!isfinite(angle)) return 0;
    const double half=angle/2;
    /*
     * 弧长转换成弦长：chord=distance*sin(angle/2)/(angle/2)。
     * 转角近零时使用 sinc 的极限展开，避免 0/0。
     */
    const double sinc=fabs(half)<1e-6 ? 1-half*half/6 : sin(half)/half;
    const double chord=distance*sinc;
    const double middle=pose->yaw+half;
    Pose2D next={pose->x+chord*cos(middle),pose->y+chord*sin(middle),
                 remainder(pose->yaw+angle,2*SENSOR_PI)};
    if (!isfinite(next.x) || !isfinite(next.y) || !isfinite(next.yaw)) return 0;
    *pose=next;
    return 1;
}
int main(void)
{
    Pose2D pose={0,0,0};
    int ok=odometry_update(&pose,1,1,0.5);
    assert(ok && fabs(pose.x-1)<1e-12 && pose.y==0);
    ok=odometry_update(&pose,-1,-1,0.5);
    assert(ok && fabs(pose.x)<1e-12); /* 倒车回到原点。 */

    pose=(Pose2D){0,0,0};
    ok=odometry_update(&pose,-SENSOR_PI/8,SENSOR_PI/8,0.5);
    assert(ok && pose.x==0 && pose.y==0 && fabs(pose.yaw-SENSOR_PI/2)<1e-12);

    /* 车体中心沿半径 1m 左转 90°：内轮半径 .75，外轮半径 1.25。 */
    pose=(Pose2D){0,0,0};
    ok=odometry_update(&pose,0.75*SENSOR_PI/2,1.25*SENSOR_PI/2,0.5);
    assert(ok && fabs(pose.x-1)<1e-12 && fabs(pose.y-1)<1e-12);
    printf("quarter circle: x=%.3f y=%.3f yaw=%.3f deg\n",pose.x,pose.y,pose.yaw*180/SENSOR_PI);
    Pose2D saved=pose;
    ok=odometry_update(&pose,1,1,0);
    assert(!ok && pose.x==saved.x);
    puts("differential_odometry: PASS");
    return 0;
}

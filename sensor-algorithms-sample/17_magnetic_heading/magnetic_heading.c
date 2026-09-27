/*
 * 磁力计倾斜补偿：先用 roll/pitch 消除倾斜，再求磁航向。
 * 右手 X前/Y左/Z上；yaw 正方向朝左（世界北-西-上约定）。
 * 输出磁北参考的 yaw，不是罗盘顺时针角，也没有加入磁偏角。
 */
#include "../common/math3d.h"
#include <assert.h>
#include <stdio.h>

static int magnetic_heading(Vec3 calibrated_mag, double roll, double pitch, double *yaw)
{
    if (yaw==NULL || !vec_finite(calibrated_mag) ||
        !isfinite(roll) || !isfinite(pitch) || fabs(pitch)>SENSOR_PI/2 ||
        fabs(cos(pitch))<1e-6) return 0;
    Vec3 normalized;
    if (!vec_unit(calibrated_mag,&normalized)) return 0;
    /*
     * R_y(pitch)*R_x(roll) 只把倾斜扶正，不额外施加未知 yaw。
     * 扶正后水平磁场在机头左侧为 +Y，对应机头偏向磁北的右侧。
     */
    Vec3 horizontal=quat_rotate(quat_from_euler(roll,pitch,0),normalized);
    if (hypot(horizontal.x,horizontal.y)<1e-6) return 0;
    *yaw=atan2(-horizontal.y,horizontal.x);
    return 1;
}
int main(void)
{
    double heading=0;
    /* 合成有倾角、有磁倾角的测量，检验补偿后能恢复已知航向。 */
    const double roll=0.3,pitch=-0.2,yaw=0.8;
    Quat truth=quat_from_euler(roll,pitch,yaw);
    Vec3 field_world={0.8,0,-0.6};
    Vec3 field_body=quat_rotate(quat_conjugate(truth),field_world);
    int ok=magnetic_heading(field_body,roll,pitch,&heading);
    assert(ok && fabs(heading-yaw)<1e-12);
    printf("magnetic yaw=%.3f deg\n",heading*180/SENSOR_PI);
    ok=magnetic_heading((Vec3){0,0,1},0,0,&heading);
    assert(!ok); /* 没有水平磁场分量时航向不确定。 */
    ok=magnetic_heading(field_body,0,SENSOR_PI/2,&heading);
    assert(!ok);
    puts("magnetic_heading: PASS");
    return 0;
}

/*
 * Pure Pursuit几何跟踪：最近路径投影+向前弧长前视点+圆弧曲率。
 * 差速/独轮车模型，输出线速度m/s与角速度rad/s。
 * 包含目标落在车后时的原地转向恢复；不直接输出阿克曼转向角。
 */
#include "../common/grid_astar.h"
#include "../common/pid_core.h"
#include <assert.h>
#include <stdio.h>
#define PATH_PI 3.14159265358979323846

typedef struct { double x,y; } Point;
typedef struct { double x,y,yaw; } RobotPose;
typedef struct { double progress; } PursuitState; /* 沿折线累计前进的弧长。 */
typedef struct { double velocity,angular_velocity,target_x,target_y; } PursuitCommand;

/* 输入为非自交、有序、相邻点不重合的折线；配置单位m、m/s、rad/s。 */
static int pursuit_update(PursuitState *state, const Point *path, int count,
                           RobotPose pose, double lookahead, double max_speed,
                           double max_angular_speed, double tolerance, PursuitCommand *out)
{
    if (state==NULL || path==NULL || out==NULL || count<2 || count>GRID_CELLS_MAX ||
        !isfinite(pose.x) || !isfinite(pose.y) || !isfinite(pose.yaw) ||
        !isfinite(state->progress) || state->progress<0 ||
        !isfinite(lookahead) || lookahead<=0 || !isfinite(max_speed) || max_speed<=0 ||
        !isfinite(max_angular_speed) || max_angular_speed<=0 ||
        !isfinite(tolerance) || tolerance<=0) return 0;
    double lengths[GRID_CELLS_MAX]={0};
    for (int i=0;i<count;++i)
    {
        if (!isfinite(path[i].x) || !isfinite(path[i].y)) return 0;
        if (i>0)
        {
            const double segment=hypot(path[i].x-path[i-1].x,path[i].y-path[i-1].y);
            if (!isfinite(segment) || segment<1e-9) return 0;
            lengths[i]=lengths[i-1]+segment;
        }
    }
    const double total=lengths[count-1];
    if (!isfinite(total) || state->progress>total+1e-9) return 0;
    double progress=state->progress,best_distance=INFINITY;
    /* 只在未走过部分找最近投影，防止正常跟踪时索引来回跳。 */
    for (int i=0;i<count-1;++i)
    {
        if (lengths[i+1]<state->progress) continue;
        const double dx=path[i+1].x-path[i].x,dy=path[i+1].y-path[i].y;
        const double segment=lengths[i+1]-lengths[i];
        const double lower=control_clamp((state->progress-lengths[i])/segment,0,1);
        double t=((pose.x-path[i].x)*dx+(pose.y-path[i].y)*dy)/(segment*segment);
        t=control_clamp(t,lower,1);
        const double distance=hypot(path[i].x+t*dx-pose.x,path[i].y+t*dy-pose.y);
        if (distance<best_distance)
        {
            best_distance=distance; progress=lengths[i]+t*segment;
        }
    }
    const double target_progress=fmin(total,progress+lookahead);
    Point target=path[count-1];
    for (int i=0;i<count-1;++i)
        if (target_progress<=lengths[i+1])
        {
            const double t=(target_progress-lengths[i])/(lengths[i+1]-lengths[i]);
            target=(Point){path[i].x+t*(path[i+1].x-path[i].x),
                           path[i].y+t*(path[i+1].y-path[i].y)};
            break;
        }
    PursuitCommand next={0,0,target.x,target.y};
    const double goal_distance=hypot(path[count-1].x-pose.x,path[count-1].y-pose.y);
    if (goal_distance>tolerance)
    {
        const double dx=target.x-pose.x,dy=target.y-pose.y;
        /* 世界目标转到机器人局部：X前、Y左。 */
        const double local_x=cos(pose.yaw)*dx+sin(pose.yaw)*dy;
        const double local_y=-sin(pose.yaw)*dx+cos(pose.yaw)*dy;
        const double distance_squared=local_x*local_x+local_y*local_y;
        if (!isfinite(distance_squared) || distance_squared<1e-12) return 0;
        if (local_x<=0)
        {
            /* 纯追踪要求目标在前方；本差速示例先原地转向恢复。 */
            next.angular_velocity=control_clamp(2*atan2(local_y,local_x),
                                                  -max_angular_speed,max_angular_speed);
        }
        else
        {
            const double curvature=2*local_y/distance_squared;
            const double stopping_gain=1.0; /* 单位1/s，将终点距离换成速度目标。 */
            next.velocity=fmin(max_speed,stopping_gain*goal_distance);
            if (fabs(curvature)>1e-12)
                next.velocity=fmin(next.velocity,max_angular_speed/fabs(curvature));
            next.angular_velocity=next.velocity*curvature;
        }
    }
    if (!isfinite(next.velocity) || !isfinite(next.angular_velocity)) return 0;
    state->progress=progress; *out=next;
    return 1;
}
static void robot_step(RobotPose *pose, PursuitCommand input, double dt)
{
    const double next_yaw=pose->yaw+input.angular_velocity*dt;
    if (fabs(input.angular_velocity)<1e-9)
    {
        pose->x+=input.velocity*cos(pose->yaw)*dt;
        pose->y+=input.velocity*sin(pose->yaw)*dt;
    }
    else
    {
        const double radius=input.velocity/input.angular_velocity;
        pose->x+=radius*(sin(next_yaw)-sin(pose->yaw));
        pose->y+=radius*(cos(pose->yaw)-cos(next_yaw));
    }
    pose->yaw=remainder(next_yaw,2*PATH_PI);
}
int main(void)
{
    const Point diagonal[]={{0,0},{2,2}};
    PursuitState state={0};
    PursuitCommand command;
    int ok=pursuit_update(&state,diagonal,2,(RobotPose){0,0,0},sqrt(2.0),0.6,1.5,0.05,&command);
    assert(ok && fabs(command.velocity-0.6)<1e-12 && fabs(command.angular_velocity-0.6)<1e-12);
    state=(PursuitState){0};
    ok=pursuit_update(&state,diagonal,2,(RobotPose){0,0,PATH_PI},1,0.6,1.5,0.05,&command);
    assert(ok && command.velocity==0 && fabs(command.angular_velocity)<=1.5);
    double saved=state.progress;
    ok=pursuit_update(&state,diagonal,2,(RobotPose){0,0,0},0,0.6,1.5,0.05,&command);
    assert(!ok && state.progress==saved);

    /* A*输出每格中心 -> 折线 -> Pure Pursuit；每格为1m。 */
    Grid grid={8,8,{0}};
    for (int row=0;row<6;++row) grid.blocked[row*grid.cols+3]=1;
    GridPath cells;
    ok=astar_find(&grid,1*grid.cols+1,6*grid.cols+6,&cells);
    assert(ok==1);
    Point path[GRID_CELLS_MAX];
    for (int i=0;i<cells.length;++i)
        path[i]=(Point){cells.cells[i]%grid.cols+0.5,cells.cells[i]/grid.cols+0.5};
    RobotPose robot={path[0].x,path[0].y,0};
    state=(PursuitState){0};
    int reached=0,steps=0;
    for (;steps<3000;++steps)
    {
        ok=pursuit_update(&state,path,cells.length,robot,0.6,0.6,1.5,0.05,&command);
        assert(ok && command.velocity>=0 && command.velocity<=0.6 &&
               fabs(command.angular_velocity)<=1.5+1e-12);
        robot_step(&robot,command,0.02);
        /* 只检验本示例的点机器人；真实车体需要障碍膨胀和连续碰撞检查。 */
        int row=(int)floor(robot.y),col=(int)floor(robot.x);
        assert(row>=0 && row<grid.rows && col>=0 && col<grid.cols);
        assert(!grid.blocked[row*grid.cols+col]);
        if (hypot(robot.x-path[cells.length-1].x,robot.y-path[cells.length-1].y)<=0.05)
        {
            reached=1; break;
        }
    }
    assert(reached);
    ok=pursuit_update(&state,path,cells.length,robot,0.6,0.6,1.5,0.05,&command);
    assert(ok && command.velocity==0 && command.angular_velocity==0);
    printf("A*+Pure Pursuit: x=%.5f y=%.5f steps=%d\n",robot.x,robot.y,steps+1);
    puts("pure_pursuit: PASS");
    return 0;
}

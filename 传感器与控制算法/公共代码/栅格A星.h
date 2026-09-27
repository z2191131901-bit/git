/*
 * 小型四邻接栅格A*：每步代价1，Manhattan启发式。
 * 固定上限16x16，不动态申请内存；线性扫描open集合，突出逻辑可读性。
 * row代表世界Y正向，col代表X正向；blocked非0表示障碍。
 */
#ifndef SENSOR_GRID_ASTAR_H
#define SENSOR_GRID_ASTAR_H
#include <stddef.h>
#include <stdlib.h>
#define GRID_SIDE_MAX 16
#define GRID_CELLS_MAX (GRID_SIDE_MAX*GRID_SIDE_MAX)

typedef struct { int rows,cols; unsigned char blocked[GRID_CELLS_MAX]; } Grid;
typedef struct { int length,cells[GRID_CELLS_MAX]; } GridPath;

static inline int grid_valid(const Grid *grid)
{
    return grid!=NULL && grid->rows>0 && grid->cols>0 &&
        grid->rows<=GRID_SIDE_MAX && grid->cols<=GRID_SIDE_MAX;
}
static inline int grid_heuristic(int cell, int goal, int cols)
{
    return abs(cell/cols-goal/cols)+abs(cell%cols-goal%cols);
}
/*
 * start/goal为row*cols+col的线性下标。
 * 返回1成功，0无路径/端点受阻，-1参数无效。失败时不改output。
 */
static inline int astar_find(const Grid *grid, int start, int goal, GridPath *output)
{
    if (!grid_valid(grid) || output==NULL) return -1;
    const int count=grid->rows*grid->cols;
    if (start<0 || start>=count || goal<0 || goal>=count) return -1;
    if (grid->blocked[start] || grid->blocked[goal]) return 0;
    int distance[GRID_CELLS_MAX],parent[GRID_CELLS_MAX];
    unsigned char closed[GRID_CELLS_MAX]={0};
    for (int i=0;i<count;++i) { distance[i]=1000000; parent[i]=-1; }
    distance[start]=0;
    const int dr[4]={0,1,0,-1},dc[4]={1,0,-1,0};
    for (int expansion=0;expansion<count;++expansion)
    {
        int current=-1,best=1000000;
        for (int cell=0;cell<count;++cell)
        {
            if (closed[cell] || distance[cell]==1000000) continue;
            const int score=distance[cell]+grid_heuristic(cell,goal,grid->cols);
            if (score<best) { best=score; current=cell; }
        }
        if (current<0) return 0;
        if (current==goal)
        {
            /* 从终点沿parent反向找回起点，再反转为执行顺序。 */
            GridPath next={0,{0}};
            for (int cell=goal;cell!=-1;cell=parent[cell])
            {
                if (next.length>=count) return -1;
                next.cells[next.length++]=cell;
            }
            for (int i=0;i<next.length/2;++i)
            {
                int temp=next.cells[i];
                next.cells[i]=next.cells[next.length-1-i];
                next.cells[next.length-1-i]=temp;
            }
            *output=next;
            return 1;
        }
        closed[current]=1;
        const int row=current/grid->cols,col=current%grid->cols;
        for (int direction=0;direction<4;++direction)
        {
            const int r=row+dr[direction],c=col+dc[direction];
            if (r<0 || r>=grid->rows || c<0 || c>=grid->cols) continue;
            const int neighbour=r*grid->cols+c;
            if (closed[neighbour] || grid->blocked[neighbour]) continue;
            const int candidate=distance[current]+1;
            if (candidate<distance[neighbour])
            {
                distance[neighbour]=candidate;
                parent[neighbour]=current;
            }
        }
    }
    return 0;
}
#endif

/* A*示例与独立BFS最短步数对照，实际实现见公共代码/栅格A星.h。 */
#include "../公共代码/栅格A星.h"
#include <assert.h>
#include <stdio.h>

/* 单位边权图上BFS可作独立最短路径基准；仅供本示例合法网格使用。 */
static int bfs_distance(const Grid *grid, int start, int goal)
{
    int queue[GRID_CELLS_MAX],distance[GRID_CELLS_MAX],head=0,tail=0;
    for (int i=0;i<grid->rows*grid->cols;++i) distance[i]=-1;
    queue[tail++]=start; distance[start]=0;
    const int dr[4]={-1,1,0,0},dc[4]={0,0,-1,1};
    while (head<tail)
    {
        int cell=queue[head++];
        if (cell==goal) return distance[cell];
        for (int d=0;d<4;++d)
        {
            int r=cell/grid->cols+dr[d],c=cell%grid->cols+dc[d];
            if (r<0 || r>=grid->rows || c<0 || c>=grid->cols) continue;
            int next=r*grid->cols+c;
            if (!grid->blocked[next] && distance[next]<0)
            {
                distance[next]=distance[cell]+1; queue[tail++]=next;
            }
        }
    }
    return -1;
}
int main(void)
{
    Grid grid={8,8,{0}};
    /* 竖墙只有最高row=7处能绕过。 */
    for (int row=0;row<7;++row) grid.blocked[row*grid.cols+3]=1;
    int start=1*grid.cols+1,goal=1*grid.cols+6;
    GridPath path;
    int ok=astar_find(&grid,start,goal,&path);
    assert(ok==1 && path.cells[0]==start && path.cells[path.length-1]==goal);
    int oracle=bfs_distance(&grid,start,goal);
    assert(oracle==17 && path.length-1==oracle);
    for (int i=0;i<path.length;++i)
    {
        assert(!grid.blocked[path.cells[i]]);
        if (i>0) assert(grid_heuristic(path.cells[i],path.cells[i-1],grid.cols)==1);
    }
    printf("A*: steps=%d BFS=%d\n",path.length-1,oracle);
    grid.blocked[7*grid.cols+3]=1;
    int saved_length=path.length;
    ok=astar_find(&grid,start,goal,&path);
    assert(ok==0 && path.length==saved_length && bfs_distance(&grid,start,goal)==-1);
    ok=astar_find(&grid,start,start,&path);
    assert(ok==1 && path.length==1);
    ok=astar_find(&grid,-1,goal,&path);
    assert(ok==-1);
    grid.blocked[start]=1;
    ok=astar_find(&grid,start,goal,&path);
    assert(ok==0);
    puts("astar: PASS");
    return 0;
}

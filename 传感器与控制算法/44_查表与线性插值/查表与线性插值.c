/* 非等间隔表：X严格递增，Y可以不单调；区间外使用端点值。 */
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
static int interpolate(const double *x,const double *y,size_t n,
                       double query,double *output)
{
    if(!x || !y || !output || n<2 || !isfinite(query)) return 0;
    /* 每次检查整张表，便于教学；嵌入式中可在初始化时检查一次。 */
    for(size_t i=0;i<n;i++)
        if(!isfinite(x[i]) || !isfinite(y[i]) ||
           (i>0 && x[i]<=x[i-1])) return 0;
    if(query<=x[0]) { *output=y[0]; return 1; }
    if(query>=x[n-1]) { *output=y[n-1]; return 1; }
    size_t low=0,high=n-1;
    /* 二分查找包含 query 的相邻两个节点，不是寻找最近的单个节点。 */
    while(high-low>1) {
        size_t middle=low+(high-low)/2;
        if(query<x[middle]) high=middle;
        else low=middle;
    }
    double span=x[high]-x[low];
    if(!isfinite(span)) return 0;
    double t=(query-x[low])/span;
    double result=(1-t)*y[low]+t*y[high];
    if(!isfinite(result)) return 0;
    *output=result;
    return 1;
}
int main(void)
{
    const double x[]={0,10,30},y[]={0,100,200},bad[]={0,10,10};
    double output=0;
    assert(interpolate(x,y,3,20,&output) && output==150);
    assert(interpolate(x,y,3,10,&output) && output==100);
    assert(interpolate(x,y,3,-1,&output) && output==0);
    assert(interpolate(x,y,3,99,&output) && output==200);
    assert(!interpolate(bad,y,3,5,&output) && output==200);
    assert(!interpolate(x,y,1,5,&output));
    puts("interpolation: PASS, x=20 -> y=150");
    return 0;
}

/* 限幅滤波：单次变化超过门限就保留上次有效值，而不是裁剪输入幅值。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
typedef struct { double value; int initialized; } Filter;
/* 1=接受，0=跳变拒绝，-1=参数错误；后两种情况不修改状态。 */
static int update(Filter *f,double input,double max_delta)
{
    if (!f || !isfinite(input) || !isfinite(max_delta) || max_delta<0)
        return -1;
    if (f->initialized && fabs(input-f->value)>max_delta) return 0;
    f->value=input; f->initialized=1;
    return 1;
}
int main(void)
{
    Filter f={0,0};
    assert(update(&f,10,2)==1);
    assert(update(&f,11,2)==1);
    assert(update(&f,100,2)==0 && f.value==11);
    assert(update(&f,13,2)==1); /* 等于门限允许通过。 */
    assert(update(&f,NAN,2)==-1 && f.value==13);
    /* 真值永久跳到 20，也会一直被拒绝：需要上层决定重新初始化。 */
    for(int i=0;i<5;i++) assert(update(&f,20,2)==0);
    f.initialized=0;
    assert(update(&f,20,2)==1);
    puts("limit filter: PASS, spike rejected; explicit reset accepts new level");
    return 0;
}

/* 固定周期采样：连续 N 次与当前稳定状态不同，才确认状态改变。
 * 输入统一成 1=按下，0=松开；低电平有效的引脚要先取反。 */
#include <assert.h>
#include <stdio.h>
typedef struct { int stable; unsigned count; } Button;
/* 返回 1=按下事件，-1=松开事件，0=无事件，2=参数错误。 */
static int sample(Button *b,int raw,unsigned required)
{
    if(!b || (raw!=0 && raw!=1) || required==0 ||
       (b->stable!=0 && b->stable!=1) || b->count>=required) return 2;
    if(raw==b->stable) { b->count=0; return 0; }
    /* 二值输入中，与 stable 不同只可能是另一个值，无需候选变量。 */
    if(++b->count<required) return 0;
    b->stable=raw; b->count=0;
    return raw ? 1 : -1;
}
int main(void)
{
    Button b={0,0};
    const int raw[]={1,0,1,1,0,1,1,1,1,0,1,0,0,0};
    int presses=0,releases=0;
    for(int i=0;i<14;i++) {
        int event=sample(&b,raw[i],3);
        presses+=(event==1); releases+=(event==-1);
    }
    assert(presses==1 && releases==1 && b.stable==0);
    assert(sample(&b,2,3)==2 && b.count==0);
    assert(sample(&b,1,0)==2);
    assert(sample(&b,1,1)==1);
    puts("debounce: PASS, one press and one release");
    return 0;
}

/* 双阈值开关：高于上阈值置1，低于下阈值置0，中间保持历史。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
static int hysteresis(int *on,double input,double low,double high)
{
    if(!on || (*on!=0 && *on!=1) || !isfinite(input) ||
       !isfinite(low) || !isfinite(high) || low>=high) return 0;
    if(input>=high) *on=1;
    else if(input<=low) *on=0;
    return 1;
}
int main(void)
{
    int on=0, changes=0;
    const double samples[]={29,31,29.5,30.5,31,27};
    const int expected[]={0,1,1,1,1,0};
    for(int i=0;i<6;i++) {
        int old=on;
        assert(hysteresis(&on,samples[i],28,31));
        assert(on==expected[i]);
        changes+=(old!=on);
    }
    assert(changes==2);
    assert(!hysteresis(&on,NAN,28,31) && on==0);
    assert(!hysteresis(&on,30,31,28));
    printf("hysteresis: PASS, transitions=%d\n",changes);
    return 0;
}

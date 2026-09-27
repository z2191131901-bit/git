/* 连续缩放死区：摇杆输入[-1,1]，死区内输出0，区外重新映射到[-1,1]。 */
#include <assert.h>
#include <math.h>
#include <stdio.h>
static int deadzone(double input,double width,double *output)
{
    if(!output || !isfinite(input) || !isfinite(width) ||
       width<0 || width>=1 || fabs(input)>1) return 0;
    const double magnitude=fabs(input);
    *output=magnitude<=width ? 0 : copysign((magnitude-width)/(1-width),input);
    return 1;
}
int main(void)
{
    double y=99;
    assert(deadzone(.1,.1,&y) && y==0);
    assert(deadzone(.55,.1,&y) && fabs(y-.5)<1e-12);
    assert(deadzone(-.55,.1,&y) && fabs(y+.5)<1e-12);
    assert(deadzone(1,.1,&y) && y==1);
    assert(deadzone(-1,.1,&y) && y==-1);
    assert(!deadzone(2,.1,&y) && y==-1);
    assert(!deadzone(0,1,&y));
    printf("deadzone: PASS, input 0.55 -> output 0.50\n");
    return 0;
}

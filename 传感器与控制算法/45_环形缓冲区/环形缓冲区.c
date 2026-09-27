/* 固定容量FIFO：满时拒绝新数据，不覆盖旧数据；单线程教学版本。 */
#include <assert.h>
#include <stddef.h>
#include <stdio.h>
enum { CAPACITY=4 };
typedef struct { int data[CAPACITY]; size_t head,tail,count; } Ring;
static int push(Ring *r,int value)
{
    if(!r || r->count==CAPACITY) return 0;
    r->data[r->head]=value;
    r->head=(r->head+1)%CAPACITY; /* 写下一个位置，到末尾后回到0。 */
    r->count++;
    return 1;
}
static int pop(Ring *r,int *value)
{
    if(!r || !value || r->count==0) return 0;
    *value=r->data[r->tail];
    r->tail=(r->tail+1)%CAPACITY;
    r->count--;
    return 1;
}
int main(void)
{
    Ring r={{0},0,0,0};
    int value=-1;
    assert(!pop(&r,&value) && value==-1);
    for(int i=1;i<=4;i++) assert(push(&r,i));
    assert(!push(&r,99) && r.count==4);
    assert(pop(&r,&value) && value==1);
    assert(pop(&r,&value) && value==2);
    assert(push(&r,5) && push(&r,6)); /* 写指针回绕后仍保持FIFO顺序。 */
    for(int i=3;i<=6;i++) assert(pop(&r,&value) && value==i);
    assert(r.count==0 && !pop(&r,&value));
    puts("ring buffer: PASS, wraparound preserves FIFO");
    return 0;
}

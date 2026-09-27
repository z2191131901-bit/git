/* 教学电机流程：待机 -> 启动中 -> 运行 -> 待机；故障优先。
 * 每次tick输入当前条件，输出enable代表是否允许电机驱动。 */
#include <assert.h>
#include <stdio.h>
typedef enum { IDLE,STARTING,RUNNING,FAULT } State;
typedef struct { State state; unsigned elapsed; int enable; } Machine;
typedef struct { int start,stop,ready,fault,reset; } Input;
static void tick(Machine *m,Input in)
{
    const unsigned timeout_ticks=3;
    /* 无论处于什么状态，当前故障都优先于启动、停止和复位。 */
    if(in.fault) { m->state=FAULT; m->elapsed=0; }
    else switch(m->state) {
    case IDLE:
        if(in.start && !in.stop) { m->state=STARTING; m->elapsed=0; }
        break;
    case STARTING:
        if(in.stop) { m->state=IDLE; m->elapsed=0; }
        else if(in.ready) { m->state=RUNNING; m->elapsed=0; }
        else if(++m->elapsed>=timeout_ticks) m->state=FAULT;
        break;
    case RUNNING:
        if(in.stop) m->state=IDLE;
        break;
    case FAULT:
        /* 故障已消失且操作者复位才回待机；这一拍不会直接重启。 */
        if(in.reset) { m->state=IDLE; m->elapsed=0; }
        break;
    default: m->state=FAULT; m->elapsed=0; break;
    }
    m->enable=(m->state==STARTING || m->state==RUNNING);
}
int main(void)
{
    Machine m={IDLE,0,0};
    tick(&m,(Input){1,0,0,0,0}); assert(m.state==STARTING && m.enable);
    tick(&m,(Input){0,0,1,0,0}); assert(m.state==RUNNING);
    tick(&m,(Input){0,0,0,1,1}); assert(m.state==FAULT && !m.enable);
    tick(&m,(Input){1,0,0,0,0}); assert(m.state==FAULT);
    tick(&m,(Input){1,0,0,0,1}); assert(m.state==IDLE && !m.enable);
    tick(&m,(Input){1,1,0,0,0}); assert(m.state==IDLE);
    tick(&m,(Input){1,0,0,0,0});
    for(int i=0;i<3;i++) tick(&m,(Input){0,0,0,0,0});
    assert(m.state==FAULT && !m.enable);
    tick(&m,(Input){0,0,0,0,1});
    tick(&m,(Input){1,0,0,0,0});
    tick(&m,(Input){0,1,1,0,0}); assert(m.state==IDLE);
    puts("state machine: PASS, fault priority, reset, stop and timeout");
    return 0;
}

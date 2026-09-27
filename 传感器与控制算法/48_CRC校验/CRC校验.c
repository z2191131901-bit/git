/* CRC-16/MODBUS：初始FFFF，反射多项式A001，最终不异或。
 * 算法数值和发送字节顺序分开：Modbus RTU先发低字节。 */
#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
static uint16_t crc_update(uint16_t crc,const uint8_t *data,size_t size)
{
    /* 调用约定：size>0时，data必须指向至少size字节的有效内存。 */
    for(size_t i=0;i<size;i++) {
        crc^=data[i];
        for(unsigned bit=0;bit<8;bit++) {
            if(crc&1u) crc=(uint16_t)((crc>>1)^0xA001u);
            else crc=(uint16_t)(crc>>1);
        }
    }
    return crc;
}
int main(void)
{
    const uint8_t text[]="123456789";
    uint16_t crc=crc_update(0xFFFFu,text,9);
    assert(crc==0x4B37u); /* 公开标准校验向量，不由本函数生成期望值。 */
    assert(crc_update(0xFFFFu,NULL,0)==0xFFFFu);
    uint16_t partial=crc_update(0xFFFFu,text,4);
    assert(crc_update(partial,text+4,5)==crc);
    uint8_t frame[]={1,3,0,0,0,10,0xC5,0xCD};
    assert(crc_update(0xFFFFu,frame,6)==0xCDC5u);
    assert(crc_update(0xFFFFu,frame,8)==0); /* 包括CRC字段后余数为0。 */
    frame[3]^=1;
    assert(crc_update(0xFFFFu,frame,8)!=0);
    printf("CRC-16/MODBUS: PASS, check=0x%04X\n",(unsigned)crc);
    return 0;
}

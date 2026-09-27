/*
 * 公共数学基础：右手系，Hamilton 乘法，q=[w,x,y,z] 表示机体到世界。
 * 世界 Z 向上；静止水平加速度计比力为 [0,0,+g]。
 * 角度使用 rad。为便于核对公式，本批用 double 中间量和状态。
 * static inline 允许头文件被各个独立示例直接包含，不需要链接额外 .c。
 */
#ifndef SENSOR_MATH3D_H
#define SENSOR_MATH3D_H
#include <math.h>
#include <stddef.h>
#define SENSOR_PI 3.14159265358979323846

typedef struct { double x, y, z; } Vec3;
typedef struct { double w, x, y, z; } Quat;

static inline int vec_finite(Vec3 v)
{
    return isfinite(v.x) && isfinite(v.y) && isfinite(v.z);
}
static inline double vec_norm(Vec3 v)
{
    return hypot(hypot(v.x, v.y), v.z);
}
static inline Vec3 vec_scale(Vec3 v, double k)
{
    Vec3 result = {v.x*k, v.y*k, v.z*k};
    return result;
}
static inline Vec3 vec_add(Vec3 a, Vec3 b)
{
    Vec3 result = {a.x+b.x, a.y+b.y, a.z+b.z};
    return result;
}
static inline Vec3 vec_cross(Vec3 a, Vec3 b)
{
    /* 叉积的顺序重要：a×b 与 b×a 方向相反。 */
    Vec3 result = {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
    return result;
}
static inline int vec_unit(Vec3 input, Vec3 *output)
{
    const double length = vec_norm(input);
    if (output == NULL || !vec_finite(input) || !isfinite(length) || length < 1e-12)
        return 0;
    *output = vec_scale(input, 1.0/length);
    return 1;
}
static inline int quat_normalize(Quat *q)
{
    if (q == NULL) return 0;
    const double length = hypot(hypot(q->w,q->x),hypot(q->y,q->z));
    if (!isfinite(length) || length < 1e-12) return 0;
    q->w /= length; q->x /= length; q->y /= length; q->z /= length;
    return 1;
}
static inline Quat quat_conjugate(Quat q)
{
    /* 只有单位四元数的逆才等于共轭。 */
    Quat result = {q.w,-q.x,-q.y,-q.z};
    return result;
}
static inline Quat quat_multiply(Quat a, Quat b)
{
    Quat result = {
        a.w*b.w-a.x*b.x-a.y*b.y-a.z*b.z,
        a.w*b.x+a.x*b.w+a.y*b.z-a.z*b.y,
        a.w*b.y-a.x*b.z+a.y*b.w+a.z*b.x,
        a.w*b.z+a.x*b.y-a.y*b.x+a.z*b.w
    };
    return result;
}
static inline Vec3 quat_rotate(Quat unit_q, Vec3 v)
{
    /* v 当作 [0,v]，q*[0,v]*共轭(q) 的向量部分就是旋转后的向量。
     * 调用者必须先保证 q 为单位四元数；此处不偷偷修改输入。 */
    Quat pure = {0,v.x,v.y,v.z};
    Quat result = quat_multiply(quat_multiply(unit_q,pure),quat_conjugate(unit_q));
    Vec3 rotated = {result.x,result.y,result.z};
    return rotated;
}
static inline Quat quat_from_euler(double roll, double pitch, double yaw)
{
    /* Z-Y-X：q = qz(yaw)*qy(pitch)*qx(roll)。输入必须为有限弧度。 */
    const double cr=cos(roll/2), sr=sin(roll/2);
    const double cp=cos(pitch/2), sp=sin(pitch/2);
    const double cy=cos(yaw/2), sy=sin(yaw/2);
    Quat q = {cr*cp*cy+sr*sp*sy, sr*cp*cy-cr*sp*sy,
              cr*sp*cy+sr*cp*sy, cr*cp*sy-sr*sp*cy};
    return q;
}
static inline int quat_step(Quat *q, Vec3 rate, double dt)
{
    if (q == NULL || !vec_finite(rate) || !isfinite(dt) || dt <= 0) return 0;
    Quat previous = *q;
    if (!quat_normalize(&previous)) return 0;
    const double speed = vec_norm(rate);
    const double half_angle = speed*dt/2;
    if (!isfinite(speed) || !isfinite(half_angle)) return 0;

    /* sin(theta/2)/speed 在 speed=0 的极限为 dt/2，避免除零。
     * 小角度泰勒展开也避免两个很小的数相除带来的数值问题。 */
    const double factor = fabs(half_angle) < 1e-8
        ? dt/2*(1-half_angle*half_angle/6) : sin(half_angle)/speed;
    Quat delta = {cos(half_angle),rate.x*factor,rate.y*factor,rate.z*factor};
    Quat next = quat_multiply(previous,delta); /* 机体角速度增量右乘。 */
    if (!quat_normalize(&next)) return 0;
    *q = next;
    return 1;
}
static inline double clamp_value(double x, double low, double high)
{
    return fmax(low,fmin(high,x));
}
#endif

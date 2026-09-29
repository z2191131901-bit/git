
/*
 * 运行管理状态：校准与故障锁存保留，停机时只清空控制历史。
 * gate_enable需要由硬件适配层落实；pwm为0不等于安全关断。
 */
typedef struct {
    FocZero zero;
    FocGuard guard;
    CurrentLoop current;
    PidState speed;
    double ramped_speed, iq_reference;
    int divider;
    PwmResult pwm;
} FocDrive;

static void drive_clear(FocDrive *s)
{
    s->current = (CurrentLoop){0};
    s->speed = (PidState){0};
    s->ramped_speed = s->iq_reference = 0.0;
    s->divider = 0;
    s->pwm = (PwmResult){0};
}

/* 输入raw已换算为A，angle为机械rad，speed为机械rad/s。 */
static int drive_tick(FocDrive *s, MotorAbc raw, double angle, double speed,
                      double bus, double target, int encoder_valid,
                      int enable, int reset)
{
    const FocGuardLimits limits = {8.0, 30.0, 6.0, 240.0};
    const PidConfig cfg = {0.14916666666666667, 1.6666666666666667,
                           0.0, 0.0, -4.0, 4.0};
    MotorAbc measured = raw;
    double theta = 0.0;
    int valid = encoder_valid && isfinite(target) &&
        foc_electrical_angle(angle, 4, 1, 0.0, &theta);
    if (s->zero.ready && !foc_zero_correct(&s->zero, raw, &measured))
        valid = 0;
    if (!foc_guard_update(&s->guard, limits, enable, reset, measured,
                          bus, speed, valid, s->zero.ready)) {
        drive_clear(s);
        return 0;
    }
    /* 电流环10kHz，速度环1kHz；使用上一电流拍的饱和反馈。 */
    if (s->divider == 0) {
        if (!foc_ramp(s->ramped_speed, target, 300.0, 10.0 * DT,
                       &s->ramped_speed) ||
            !foc_speed_pi(&s->speed, &cfg, s->ramped_speed, speed,
                          10.0 * DT, s->current.voltage_limited,
                          s->current.q_voltage_deficit, &s->iq_reference))
            goto failure;
    }
    MotorDq0 reference;
    if (!foc_current_circle(0.0, s->iq_reference, 4.0, &reference) ||
        !foc_step(&s->current, measured, theta, speed,
                   reference.d, reference.q, bus, 1, &s->pwm))
        goto failure;
    s->divider = (s->divider + 1) % 10;
    return 1;
failure:
    s->guard.faults |= FOC_BAD_INPUT;
    s->guard.gate_enable = 0;
    drive_clear(s);
    return 0;
}

static void check_input_and_guard(void)
{
    FocZero zero = {0};
    MotorAbc offset = {0.1, -0.05, 0.02}, corrected;
    assert(!foc_zero_sample(&zero, offset, 0, 1));
    assert(!foc_zero_sample(&zero, offset, 1, 0));
    assert(!foc_zero_sample(&zero, (MotorAbc){NAN,0,0}, 1, 1));
    assert(zero.count == 0);
    for (int k = 0; k < FOC_ZERO_SAMPLES; ++k) {
        double noise = 0.02 * cos(2.0 * MOTOR_PI * k / FOC_ZERO_SAMPLES);
        assert(foc_zero_sample(&zero,
            (MotorAbc){offset.a+noise,offset.b-noise,offset.c+noise}, 1, 1));
    }
    assert(zero.ready && fabs(zero.offset.a-offset.a)<1e-12);
    assert(foc_zero_correct(&zero, (MotorAbc){1.1,-0.55,-0.48}, &corrected));
    assert(fabs(corrected.a-1.0)<1e-12 && fabs(corrected.b+0.5)<1e-12 &&
           fabs(corrected.c+0.5)<1e-12);
    double angle = 99.0;
    assert(foc_electrical_angle(MOTOR_PI/8.0,4,1,0.0,&angle));
    assert(fabs(angle-MOTOR_PI/2.0)<1e-12);
    assert(foc_electrical_angle(MOTOR_PI/8.0,4,-1,MOTOR_PI/2.0,&angle));
    assert(fabs(angle)<1e-12);
    assert(!foc_electrical_angle(0.0,0,1,0.0,&angle) && angle==0.0);
    double ramp = 0.0;
    for (int k=0;k<1000;++k)
        assert(foc_ramp(ramp,10.0,100.0,0.001,&ramp));
    assert(fabs(ramp-10.0)<1e-12);
    assert(foc_ramp(ramp,-10.0,100.0,0.001,&ramp) && fabs(ramp-9.9)<1e-12);
    assert(!foc_ramp(ramp,NAN,100.0,0.001,&ramp) && fabs(ramp-9.9)<1e-12);
    MotorDq0 ref;
    assert(foc_current_circle(3.0,4.0,4.0,&ref));
    assert(fabs(ref.q-sqrt(7.0))<1e-12 && fabs(hypot(ref.d,ref.q)-4.0)<1e-12);
    assert(foc_current_circle(-5.0,-4.0,4.0,&ref) && ref.d==-4.0 && ref.q==0.0);

    /* 分别注入每类故障；正常测量不清除锁存，运行中reset无效。 */
    const FocGuardLimits limits = {8,30,6,240};
    const unsigned expected[] = {FOC_UNDERVOLTAGE,FOC_OVERVOLTAGE,
                                 FOC_OVERCURRENT,FOC_OVERSPEED,FOC_BAD_INPUT};
    for (int k=0;k<5;++k) {
        FocGuard guard = {0};
        MotorAbc sample = {k==2?6.1:0.0,0,0};
        double bus = k==0?7.9:k==1?30.1:24.0;
        double speed = k==3?-240.1:0.0;
        assert(!foc_guard_update(&guard,limits,1,0,sample,bus,speed,k!=4,1));
        assert(guard.faults==expected[k]);
        assert(!foc_guard_update(&guard,limits,1,1,(MotorAbc){0},24,0,1,1));
        assert(guard.faults==expected[k]);
        assert(!foc_guard_update(&guard,limits,0,1,(MotorAbc){0},24,0,1,1));
        assert(guard.faults==0 && !guard.gate_enable);
        assert(foc_guard_update(&guard,limits,1,0,(MotorAbc){0},24,0,1,1));
    }
    FocGuard guard = {0};
    assert(!foc_guard_update(&guard,limits,0,1,(MotorAbc){0},NAN,0,1,1));
    assert(guard.faults==FOC_BAD_INPUT);
    const PidConfig cfg = {0.1,1.0,0,0,-4,4};
    PidState pi = {0};
    double command;
    assert(foc_speed_pi(&pi,&cfg,1,0,0.001,1,2,&command) && pi.integral==0);
    assert(foc_speed_pi(&pi,&cfg,-1,0,0.001,1,2,&command) && pi.integral<0);
    double saved = pi.integral;
    assert(foc_speed_pi(&pi,&cfg,-1,0,0.001,1,-2,&command) && pi.integral==saved);
    assert(!foc_speed_pi(&pi,&cfg,NAN,0,0.001,0,0,&command) && pi.integral==saved);

    /* 完整链：带零偏的三相测量+机械编码器角度，运行0.8秒。 */
    FocDrive drive = {0};
    assert(!drive_tick(&drive,offset,0,0,24,60,1,1,0));
    drive.zero = zero;
    MotorState motor = {0};
    for (int k=0;k<8000;++k) {
        MotorAbc raw = measured_currents(motor);
        raw.a+=offset.a;raw.b+=offset.b;raw.c+=offset.c;
        assert(drive_tick(&drive,raw,motor.theta/POLE_PAIRS,motor.speed,
                           24,60,1,1,0));
        MotorAb0 voltage = actual_voltage(drive.pwm,24);
        motor_step(&motor,voltage.alpha,voltage.beta,0,0,4);
    }
    assert(fabs(motor.speed-60.0)<0.1);
    printf("FOC: calibrated drive final_speed=%.6f rad/s\n",motor.speed);
    /* 不模拟关断后的功率电路；这里检查控制与使能状态。 */
    assert(!drive_tick(&drive,(MotorAbc){7.1,-0.05,0.02},0,0,24,60,1,1,0));
    assert(drive.guard.faults==FOC_OVERCURRENT && !drive.guard.gate_enable);
    assert(drive.current.integral_q==0 && drive.speed.integral==0 &&
           drive.ramped_speed==0 && drive.zero.ready);
    assert(!drive_tick(&drive,offset,0,0,24,60,1,1,0));
    assert(!drive_tick(&drive,offset,0,0,24,60,1,0,1));
    assert(drive.guard.faults==0);
    assert(drive_tick(&drive,offset,0,0,24,60,1,1,0));
    assert(fabs(drive.ramped_speed-0.3)<1e-12);
    assert(!drive_tick(&drive,offset,0,0,24,60,1,0,0));
    assert(drive.speed.integral==0 && drive.ramped_speed==0);
    puts("PASS: calibration, electrical angle, ramp, current circle, latched faults and restart.");
}

/* 母线8V时80rad/s不可持续达到；恢复24V后同时把目标降为20rad/s。 */
static double check_outer_saturation(int coordinated)
{
    MotorState motor = {0};
    CurrentLoop inner = {0};
    PidState outer = {0};
    const PidConfig cfg = {0.14916666666666667,1.6666666666666667,0,0,-4,4};
    double iq = 0, iae = 0;
    int limited = 0;
    for (int k=0;k<20000;++k) {
        double target=k<10000?80:20, bus=k<10000?8:24;
        if (k%10==0) {
            if (coordinated)
                assert(foc_speed_pi(&outer,&cfg,target,motor.speed,10*DT,
                                    inner.voltage_limited,inner.q_voltage_deficit,&iq));
            else
                assert(pid_update(&outer,&cfg,target,motor.speed,0,10*DT,&iq));
        }
        PwmResult pwm;
        assert(foc_step(&inner,measured_currents(motor),motor.theta,motor.speed,
                         0,iq,bus,1,&pwm));
        limited+=pwm.limited;
        MotorAb0 v=actual_voltage(pwm,bus);
        motor_step(&motor,v.alpha,v.beta,0,0,4);
        if(k>=10000) iae+=fabs(target-motor.speed)*DT;
    }
    assert(limited>0 && fabs(motor.speed-20)<0.1);
    printf("FOC: outer coordination=%d recovery_IAE=%.9f rad final=%.6f rad/s\n",
           coordinated,iae,motor.speed);
    return iae;
}

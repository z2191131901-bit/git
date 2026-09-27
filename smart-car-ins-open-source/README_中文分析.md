# 智能车惯性导航开源代码调研与源码分析

更新日期：2026-09-24  
代码位置：本目录下 4 个独立 Git 仓库，均保留其 Git 元数据和原始许可证。

## 1. 下载结果与建议

| 项目 | 本地目录 | 主要传感器 | 核心方法 | 适合用途 | 许可证 |
|---|---|---|---|---|---|
| [MapIV/eagleye](https://github.com/MapIV/eagleye) | `eagleye/` | GNSS（位置/多普勒）、IMU、车速 | 面向地面车辆的分模块航向、速度、位置、侧滑估计 | ROS 2 实车定位、理解车辆专用导航链 | BSD 3-Clause |
| [rsasaki0109/kalman_filter_localization_ros2](https://github.com/rsasaki0109/kalman_filter_localization_ros2) | `kalman_filter_localization_ros2/` | GNSS、IMU、轮速 | 15 维误差状态 EKF，带车辆约束、延迟回放和鲁棒处理 | 最值得二次开发的 ROS 2 基线 | BSD 2-Clause |
| [jnz/INSLIB](https://github.com/jnz/INSLIB) | `INSLIB/` | GNSS、IMU、轮速、磁力计、气压计等 | 嵌入式多传感器惯性组合导航、UDU/Bierman 滤波 | MCU/嵌入式和高可靠实现参考 | AGPL-3.0 |
| [jasleon/Vehicle-State-Estimation](https://github.com/jasleon/Vehicle-State-Estimation) | `Vehicle-State-Estimation/` | IMU、GNSS、LiDAR 位姿 | 9 维 ES-EKF 教学实现 | 快速看懂预测—校正全过程 | 仓库未声明许可证，不应复制进产品 |

如果目标是“智能车上真正跑起来”，推荐以 `kalman_filter_localization_ros2` 为主线，以 Eagleye 的 GNSS 多普勒航向/速度处理作为车辆特化参考。若控制器不是 Linux/ROS，而是 MCU，则优先评估 INSLIB，但必须先接受 AGPL-3.0 的开源义务。教学仓库只用于学习公式与代码映射。

## 2. 它们共同是怎么做到的

惯性导航的基本矛盾是：IMU 更新快且不依赖外界，但积分会漂移；GNSS、LiDAR 或地图定位给出低频绝对位置，轮速给出稳定的车体前向速度。组合导航用高频 IMU 做预测，用低频外部观测不断把漂移拉回来。

### 2.1 名义状态传播

典型状态为位置、速度、姿态和 IMU 零偏：

```text
x = [p, v, q, bg, ba]
```

每个 IMU 周期先去零偏：

```text
omega = gyro - bg
a     = accel - ba
```

再完成捷联解算：角速度积分更新四元数，机体系加速度由姿态旋转到世界系并扣除重力，然后一次积分得速度、二次积分得位置：

```text
q(k+1) = q(k) ⊗ Exp(omega * dt)
v(k+1) = v(k) + (R(q) * a - g) * dt
p(k+1) = p(k) + v(k) * dt + 0.5 * (R(q) * a - g) * dt²
```

### 2.2 误差状态 EKF

工程实现通常不直接在线性空间里加减四元数，而是维护 15 维小误差：

```text
delta_x = [delta_p, delta_v, delta_theta, delta_bg, delta_ba]
```

预测阶段同时用线性化矩阵传播协方差：

```text
P- = F P F' + Q
```

收到 GNSS、轮速或 LiDAR 测量时，计算创新 `r = z - h(x)`、创新协方差 `S` 和卡尔曼增益 `K`：

```text
S = H P- H' + R
K = P- H' inv(S)
delta_x = K r
```

随后把位置、速度、零偏误差直接注入名义状态，把小角度通过指数映射注入四元数，再把误差状态复位为零。这样既保持姿态在旋转群上合法，也比直接 EKF 更适合高动态非线性系统。

### 2.3 地面车辆额外可用的信息

智能车与通用飞行器不同，可加入很强的车辆约束：

- 轮速主要观测车体 x 轴前向速度；
- 正常行驶时车体横向和竖向速度近似为 0，即非完整约束 NHC；
- 停车时速度为 0（ZUPT），角速度也应接近 0（ZIHR/ZARU），可以校准漂移；
- GNSS 天线不在 IMU 原点，必须补偿杆臂 `p_gnss = p_body + R * lever_arm`；
- GNSS 有处理与传输延迟，最好回退到测量时刻校正后再重放 IMU，而不是用“旧测量”修正“新状态”；
- GNSS 多普勒速度通常比相邻位置差分更平滑，运动时还可提供航向可观性。

## 3. 项目逐一拆解

### 3.1 Eagleye：把车辆定位拆成可观测的小问题

关键代码：

- `eagleye_core/navigation/src/velocity_estimator.cpp`：融合 GNSS 多普勒速度、RTK 位置差分速度和 IMU 累积加速度，包含异常值判断、停车判断和滑窗处理。
- `eagleye_core/navigation/src/heading.cpp`：由 ENU 速度计算 `atan2(v_east, v_north)` 得到运动航向，并与陀螺积分配合；低速时航向不可观，所以代码有速度阈值和状态管理。
- `eagleye_core/navigation/src/yaw_rate_offset*.cpp`：车辆静止或稳定运动阶段估计陀螺 z 轴偏置。
- `eagleye_core/navigation/src/trajectory.cpp`：利用速度、航向、俯仰/横滚进行轨迹递推。
- `eagleye_core/navigation/src/slip_angle.cpp`：处理车辆实际速度方向与车头方向不一致的问题。

它不是一个“大一统 EKF 文件”，而是把速度尺度、陀螺偏置、航向、侧滑、位置和高度拆成 ROS 节点流水线。优点是每个误差源都能单独诊断和调参，尤其适合汽车；缺点是模块、话题和配置多，接入成本高，对 GNSS 质量和车辆运动条件也有要求。

### 3.2 kalman_filter_localization_ros2：完整的 15 维车辆 ESKF

关键入口：

- `docs/eskf_math.md`：该项目的坐标系、状态、雅可比和测量模型规范，应先读。
- `kalman_filter_localization_core/include/kalman_filter_localization/core/ekf_estimator.hpp`：核心预测、测量更新、误差注入。
- `kalman_filter_localization_core/include/kalman_filter_localization/core/eskf_replay.hpp`：延迟测量的历史状态回退与重放。
- `kalman_filter_localization_core/include/kalman_filter_localization/core/vehicle_observability.hpp`：车辆状态可观性判断。
- `kalman_filter_localization_ros2/src/ekf_localization_component.cpp`：ROS 2 消息、TF、时间戳和参数到核心滤波器的适配层。
- `kalman_filter_localization_ros2/param/ekf.yaml`：默认噪声、门限、杆臂和功能开关。

它的名义状态是 16 个参数（四元数占 4 个），局部误差是 15 维。加速度计和陀螺仪零偏都作为随机过程在线估计。项目同时实现：

- GNSS 位置、GNSS 多普勒速度、车轮速度；
- GNSS 天线杆臂；
- NHC、ZUPT、零角速度约束；
- NIS 门限及 Huber/Cauchy 鲁棒损失，抑制跳点；
- 延迟测量 rewind/replay；
- 二阶离散化与离线 Van Loan 精确离散化；
- 固定滞后平滑、故障注入、数据集评估和单元测试。

这套代码最接近可扩展的研究/工程基线。落地时首先要统一 ENU/FLU 坐标定义、时间同步、IMU 是否含重力、传感器外参和噪声单位；这些环节任一错误，都可能比滤波算法本身造成更大的偏差。

### 3.3 INSLIB：面向嵌入式与高可靠约束

关键入口：

- `src/ins.h` / `src/ins.c`：INS 状态、输入测量、初始化、预测与各种观测融合。
- `src/nav_suite.h`：INS、AHRS、气压高度等模块的统一封装。
- `KFCore/`：UDU/Bierman-Thornton 数值稳定卡尔曼滤波核心。
- `datasets/tunnel_odometry/`：汽车驶入隧道时使用轮速桥接 GNSS 中断的示例数据。
- `tutorial/c_tutorial.md`：C/C++ 接入入口。

它强调无堆内存、无操作系统依赖、固定栈使用、C11 和单精度友好，适合 MCU。输入除 IMU/GNSS 外，还可包含磁力计、气压计、轮速、局部位置、绝对航向、ZUPT/ZARU。其工程亮点包括 GNSS 延迟补偿、杆臂、自动静止检测、传感器精度门控、测量异常判断、超过 90% 的测试覆盖和需求追踪。

与 ROS 2 ESKF 相比，它更像可嵌入产品的导航库，但 AGPL-3.0 是强 copyleft 许可证：若修改后用于网络服务，或把它链接进分发的软件，通常会触发对应源代码开放义务。商用前必须让法务确认使用方式。

### 3.4 Vehicle-State-Estimation：最短路径理解 ES-EKF

核心全部集中在 `es_ekf.py`：

1. 用 IMU 比力和角速度传播 `p/v/q`；
2. 构造 9 维误差协方差，对应位置、速度、小姿态误差；
3. GNSS 和 LiDAR 都被简化为三维位置观测；
4. LiDAR 位姿先用固定外参旋转和平移到 IMU 坐标系；
5. 测量更新得到小误差后，将位置/速度相加，将旋转小量以四元数方式注入。

它没有估计 IMU 零偏，没有杆臂、延迟补偿、异常值门控、车辆约束或实时消息系统，而且初始化直接使用真值，因此不能原样上车。当前环境实测还会在 `rotations.py` 的 `skew_symmetric()` 因新版 NumPy 对列向量的处理而报形状错误；把输入显式展平成 `(3,)` 可修复兼容性，但本仓库没有许可证，所以这里没有改动原代码。

## 4. 实车落地的数据流

```text
IMU(100~400 Hz) ──标定/轴向变换/时间同步──> ESKF 预测 ──────────────> 高频 p,v,q
                                                     ↑
GNSS(1~20 Hz) ──ENU转换/杆臂/延迟补偿/门控────────────┤
轮速(20~100 Hz) ──比例/方向/打滑判断─────────────────┤
车辆静止 ──ZUPT + ZARU───────────────────────────────┘
```

建议实施顺序：

1. 先只用 IMU + GNSS，在开阔场地验证坐标系、重力符号和时间戳；
2. 加入轮速，验证直线、转弯、倒车和打滑场景；
3. 加入 NHC、ZUPT/ZARU，并检查其触发条件不会在侧滑或坡道时误判；
4. 标定 IMU 到车体、GNSS 天线到 IMU 的外参；
5. 用 Allan 方差或静态长采样得到 IMU 白噪声与零偏随机游走，而非照抄默认参数；
6. 人工制造 GNSS 遮挡和跳点，检查协方差、NIS、重捕获和延迟回放；
7. 最后才接入 LiDAR/视觉定位，并明确其输出坐标系和真实延迟。

## 5. 验证记录与限制

- 四个仓库均在 2026-09-24 以浅克隆方式下载；INSLIB 的 `KFCore` 子模块已递归补齐。
- 下载版本：Eagleye `1095620e`；ROS 2 ESKF `0a456b24`；Vehicle-State-Estimation `3d672307`；INSLIB `28fd0a06`，KFCore `8648bad1`。
- 教学 Python 项目已尝试在当前 Windows/Python 环境运行，因上述 NumPy 维度兼容问题中止。
- 另外三个项目没有在此机器编译：Eagleye 和 ROS 2 ESKF 需要匹配的 Linux/ROS 2 环境；INSLIB 虽支持 Windows，但完整验证需要其指定编译器和依赖环境。
- “开源”不等于可以随意摘抄。BSD 项目可在保留版权和免责声明的前提下修改/分发；AGPL 项目有强开源义务；未声明许可证的仓库只能阅读，不能默认获得复制、修改和分发授权。

## 6. 最终选型

- **ROS 2 智能车研究与实车原型**：从 `kalman_filter_localization_ros2` 开始。
- **GNSS 多普勒与车辆专用定位链**：重点研究 `eagleye`。
- **裸机/MCU/嵌入式**：技术上优先研究 `INSLIB`，同时先解决许可证问题。
- **学习 ESKF 数学和最小代码**：阅读 `Vehicle-State-Estimation`，不要直接用于产品。


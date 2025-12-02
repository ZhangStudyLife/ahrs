/*******************************************************************************
 * 姿态角速度积分算法 (Angular Rate Integration)
 * 
 * 算法类型: 基础姿态传播
 * 功能说明: 通过对陀螺仪角速度积分来更新姿态四元数
 *
 * =============================================================================
 * 理论基础：
 * =============================================================================
 * 
 * 角速度积分是最基础的姿态更新方法。当刚体在空间中旋转时，其姿态的变化
 * 可以通过对角速度进行积分来获得。
 * 
 * 【四元数微分方程】
 * 四元数 q(t) 对时间的导数与角速度 ω(t) 的关系为：
 * 
 *     dq/dt = (1/2) * q ⊗ ω
 * 
 * 其中 ω 被表示为纯四元数 [0, ωx, ωy, ωz]
 * 
 * 【数值积分方法】
 * 
 * 1. 一阶欧拉法 (最简单):
 *    q(t+Δt) = q(t) + (Δt/2) * q(t) ⊗ ω(t)
 *    优点：计算量小
 *    缺点：精度低，误差累积快
 * 
 * 2. 中点法 (本实现采用):
 *    q(t+Δt) = q(t) + (Δt/2) * [q(t) + q(t+Δt)] ⊗ ω_mid
 *    优点：精度较高
 *    缺点：需要迭代或预测
 * 
 * 3. 四阶龙格-库塔法 (最精确):
 *    使用4次中间步骤计算
 *    优点：精度最高
 *    缺点：计算量大
 * 
 * 【误差来源】
 * 1. 陀螺仪零偏 (Bias): 导致姿态漂移
 * 2. 陀螺仪白噪声: 导致姿态抖动
 * 3. 数值积分截断误差: 时间步长越大误差越大
 * 
 * =============================================================================
 * 算法特点：
 * =============================================================================
 * 
 * ✓ 优点：
 *   - 实现简单，计算量最小
 *   - 短时间内精度高
 *   - 响应速度快，无延迟
 *   - 不受磁场干扰
 * 
 * ✗ 缺点：
 *   - 存在累积误差，会漂移
 *   - 需要其他传感器校正
 *   - 对陀螺仪零偏敏感
 *   - 长时间运行需要重置
 * 
 * 应用场景：
 * - 高速姿态跟踪
 * - 与其他滤波器配合使用
 * - 短时间内的姿态估计
 * 
 ******************************************************************************/

#include <math.h>

/*===========================================================================
 * 数据结构定义
 *===========================================================================*/

// 四元数结构体
typedef struct {
    float w, x, y, z;  // w是实部，xyz是虚部
} Quaternion;

// 三维向量（角速度）
typedef struct {
    float x, y, z;
} Vector3;

// Angular积分器状态
typedef struct {
    Quaternion q;       // 当前姿态四元数
    float dt;           // 采样时间间隔（秒）
    int method;         // 积分方法: 0=欧拉法, 1=中点法, 2=龙格库塔法
} AngularIntegrator;

/*===========================================================================
 * 核心数学函数
 *===========================================================================*/

/**
 * 四元数归一化
 * 
 * 由于数值误差，四元数的模长可能偏离1，需要定期归一化
 * 归一化公式：q_norm = q / ||q||
 */
void quaternion_normalize(Quaternion* q) {
    float norm = sqrtf(q->w * q->w + q->x * q->x + 
                       q->y * q->y + q->z * q->z);
    
    // 避免除零
    if (norm > 1e-9f) {
        float inv_norm = 1.0f / norm;
        q->w *= inv_norm;
        q->x *= inv_norm;
        q->y *= inv_norm;
        q->z *= inv_norm;
    } else {
        // 如果四元数接近零，重置为单位四元数
        q->w = 1.0f;
        q->x = q->y = q->z = 0.0f;
    }
}

/**
 * 四元数乘法（Hamilton乘积）
 * 
 * 四元数乘法不满足交换律！顺序很重要
 * q1 ⊗ q2 ≠ q2 ⊗ q1
 * 
 * 乘法规则：
 * (w1 + x1*i + y1*j + z1*k) ⊗ (w2 + x2*i + y2*j + z2*k)
 * 
 * 其中: i² = j² = k² = ijk = -1
 *      ij = k,  jk = i,  ki = j
 *      ji = -k, kj = -i, ik = -j
 */
Quaternion quaternion_multiply(const Quaternion* q1, const Quaternion* q2) {
    Quaternion result;
    
    // 实部 = w1*w2 - (x1*x2 + y1*y2 + z1*z2)
    result.w = q1->w * q2->w - q1->x * q2->x - 
               q1->y * q2->y - q1->z * q2->z;
    
    // 虚部 i 分量 = w1*x2 + x1*w2 + y1*z2 - z1*y2
    result.x = q1->w * q2->x + q1->x * q2->w + 
               q1->y * q2->z - q1->z * q2->y;
    
    // 虚部 j 分量 = w1*y2 - x1*z2 + y1*w2 + z1*x2
    result.y = q1->w * q2->y - q1->x * q2->z + 
               q1->y * q2->w + q1->z * q2->x;
    
    // 虚部 k 分量 = w1*z2 + x1*y2 - y1*x2 + z1*w2
    result.z = q1->w * q2->z + q1->x * q2->y - 
               q1->y * q2->x + q1->z * q2->w;
    
    return result;
}

/**
 * 将角速度向量转换为四元数导数
 * 
 * 公式：dq/dt = (1/2) * q ⊗ ω_quat
 * 
 * 其中 ω_quat = [0, ωx, ωy, ωz] 是角速度的四元数表示
 * 
 * @param q 当前四元数
 * @param omega 角速度向量（单位：rad/s）
 * @return 四元数的导数
 */
Quaternion compute_quaternion_derivative(const Quaternion* q, const Vector3* omega) {
    // 角速度四元数：[0, ωx, ωy, ωz]
    Quaternion omega_quat;
    omega_quat.w = 0.0f;
    omega_quat.x = omega->x;
    omega_quat.y = omega->y;
    omega_quat.z = omega->z;
    
    // q ⊗ ω
    Quaternion q_dot = quaternion_multiply(q, &omega_quat);
    
    // 乘以 1/2
    q_dot.w *= 0.5f;
    q_dot.x *= 0.5f;
    q_dot.y *= 0.5f;
    q_dot.z *= 0.5f;
    
    return q_dot;
}

/*===========================================================================
 * 积分算法实现
 *===========================================================================*/

/**
 * 一阶欧拉法更新四元数
 * 
 * 最简单的积分方法：
 * q(t+Δt) ≈ q(t) + Δt * dq/dt
 *         = q(t) + Δt * (1/2) * q(t) ⊗ ω(t)
 * 
 * 优点：计算量最小
 * 缺点：精度较低，误差O(Δt²)
 * 
 * 适用场景：采样频率很高（>500Hz）的系统
 */
void angular_update_euler(AngularIntegrator* integrator, const Vector3* gyro) {
    // 计算四元数导数
    Quaternion q_dot = compute_quaternion_derivative(&integrator->q, gyro);
    
    // 欧拉积分：q_new = q_old + dt * q_dot
    integrator->q.w += integrator->dt * q_dot.w;
    integrator->q.x += integrator->dt * q_dot.x;
    integrator->q.y += integrator->dt * q_dot.y;
    integrator->q.z += integrator->dt * q_dot.z;
    
    // 归一化四元数（重要！）
    quaternion_normalize(&integrator->q);
}

/**
 * 二阶中点法（梯形法）更新四元数
 * 
 * 中点法使用两个时间点的平均值：
 * q(t+Δt) ≈ q(t) + (Δt/2) * [dq/dt(t) + dq/dt(t+Δt)]
 * 
 * 实现时使用预测-校正策略：
 * 1. 预测：q_pred = q(t) + Δt * dq/dt(t)
 * 2. 校正：q(t+Δt) = q(t) + (Δt/2) * [dq/dt(t) + dq/dt_pred]
 * 
 * 精度：误差O(Δt³)，比欧拉法高一个数量级
 */
void angular_update_midpoint(AngularIntegrator* integrator, const Vector3* gyro) {
    // 1. 计算当前时刻的导数
    Quaternion q_dot_start = compute_quaternion_derivative(&integrator->q, gyro);
    
    // 2. 预测中点状态（欧拉法半步）
    Quaternion q_mid;
    float half_dt = integrator->dt * 0.5f;
    q_mid.w = integrator->q.w + half_dt * q_dot_start.w;
    q_mid.x = integrator->q.x + half_dt * q_dot_start.x;
    q_mid.y = integrator->q.y + half_dt * q_dot_start.y;
    q_mid.z = integrator->q.z + half_dt * q_dot_start.z;
    quaternion_normalize(&q_mid);
    
    // 3. 计算中点的导数
    Quaternion q_dot_mid = compute_quaternion_derivative(&q_mid, gyro);
    
    // 4. 使用中点导数更新（全步）
    integrator->q.w += integrator->dt * q_dot_mid.w;
    integrator->q.x += integrator->dt * q_dot_mid.x;
    integrator->q.y += integrator->dt * q_dot_mid.y;
    integrator->q.z += integrator->dt * q_dot_mid.z;
    
    // 5. 归一化
    quaternion_normalize(&integrator->q);
}

/**
 * 四阶龙格-库塔法（RK4）更新四元数
 * 
 * RK4是最常用的高精度数值积分方法
 * 
 * 算法步骤：
 * k1 = f(t, q)
 * k2 = f(t + Δt/2, q + Δt*k1/2)
 * k3 = f(t + Δt/2, q + Δt*k2/2)
 * k4 = f(t + Δt, q + Δt*k3)
 * q_new = q + (Δt/6) * (k1 + 2*k2 + 2*k3 + k4)
 * 
 * 精度：误差O(Δt⁵)，非常高
 * 代价：需要计算4次导数
 */
void angular_update_rk4(AngularIntegrator* integrator, const Vector3* gyro) {
    Quaternion q0 = integrator->q;
    float dt = integrator->dt;
    
    // k1 = dq/dt at t
    Quaternion k1 = compute_quaternion_derivative(&q0, gyro);
    
    // k2 = dq/dt at t + dt/2, using q + dt*k1/2
    Quaternion q_temp;
    q_temp.w = q0.w + 0.5f * dt * k1.w;
    q_temp.x = q0.x + 0.5f * dt * k1.x;
    q_temp.y = q0.y + 0.5f * dt * k1.y;
    q_temp.z = q0.z + 0.5f * dt * k1.z;
    quaternion_normalize(&q_temp);
    Quaternion k2 = compute_quaternion_derivative(&q_temp, gyro);
    
    // k3 = dq/dt at t + dt/2, using q + dt*k2/2
    q_temp.w = q0.w + 0.5f * dt * k2.w;
    q_temp.x = q0.x + 0.5f * dt * k2.x;
    q_temp.y = q0.y + 0.5f * dt * k2.y;
    q_temp.z = q0.z + 0.5f * dt * k2.z;
    quaternion_normalize(&q_temp);
    Quaternion k3 = compute_quaternion_derivative(&q_temp, gyro);
    
    // k4 = dq/dt at t + dt, using q + dt*k3
    q_temp.w = q0.w + dt * k3.w;
    q_temp.x = q0.x + dt * k3.x;
    q_temp.y = q0.y + dt * k3.y;
    q_temp.z = q0.z + dt * k3.z;
    quaternion_normalize(&q_temp);
    Quaternion k4 = compute_quaternion_derivative(&q_temp, gyro);
    
    // 组合：q_new = q + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
    integrator->q.w += (dt / 6.0f) * (k1.w + 2.0f*k2.w + 2.0f*k3.w + k4.w);
    integrator->q.x += (dt / 6.0f) * (k1.x + 2.0f*k2.x + 2.0f*k3.x + k4.x);
    integrator->q.y += (dt / 6.0f) * (k1.y + 2.0f*k2.y + 2.0f*k3.y + k4.y);
    integrator->q.z += (dt / 6.0f) * (k1.z + 2.0f*k2.z + 2.0f*k3.z + k4.z);
    
    // 归一化
    quaternion_normalize(&integrator->q);
}

/*===========================================================================
 * 公共接口函数
 *===========================================================================*/

/**
 * 初始化角速度积分器
 * 
 * @param integrator 积分器结构体指针
 * @param dt 采样时间间隔（秒），例如100Hz采样时dt=0.01
 * @param method 积分方法: 0=欧拉法, 1=中点法, 2=RK4
 * 
 * 初始姿态设为单位四元数（无旋转）
 */
void angular_init(AngularIntegrator* integrator, float dt, int method) {
    // 初始化为单位四元数（表示无旋转）
    integrator->q.w = 1.0f;
    integrator->q.x = 0.0f;
    integrator->q.y = 0.0f;
    integrator->q.z = 0.0f;
    
    integrator->dt = dt;
    integrator->method = method;
}

/**
 * 更新姿态（主函数）
 * 
 * @param integrator 积分器结构体指针
 * @param gx, gy, gz 三轴角速度（单位：rad/s）
 * 
 * 【注意事项】
 * 1. 角速度单位必须是rad/s，如果陀螺仪输出是deg/s需要转换
 *    rad/s = deg/s * (π/180) = deg/s * 0.017453
 * 
 * 2. 坐标系定义：
 *    - X轴：前进方向（Roll，横滚）
 *    - Y轴：右侧方向（Pitch，俯仰）
 *    - Z轴：向下方向（Yaw，偏航）
 * 
 * 3. 建议每隔100-1000次更新调用一次额外的归一化
 *    以补偿累积的数值误差
 */
void angular_update(AngularIntegrator* integrator, float gx, float gy, float gz) {
    Vector3 gyro = {gx, gy, gz};
    
    // 根据选择的方法进行更新
    switch (integrator->method) {
        case 0:  // 欧拉法
            angular_update_euler(integrator, &gyro);
            break;
        case 1:  // 中点法（默认）
            angular_update_midpoint(integrator, &gyro);
            break;
        case 2:  // RK4
            angular_update_rk4(integrator, &gyro);
            break;
        default:
            angular_update_midpoint(integrator, &gyro);
    }
}

/**
 * 设置当前姿态四元数
 * 
 * 用于：
 * - 设置初始姿态
 * - 从其他传感器校正姿态
 * - 融合其他滤波器的结果
 * 
 * @param integrator 积分器结构体指针
 * @param w, x, y, z 四元数分量
 */
void angular_set_quaternion(AngularIntegrator* integrator, 
                           float w, float x, float y, float z) {
    integrator->q.w = w;
    integrator->q.x = x;
    integrator->q.y = y;
    integrator->q.z = z;
    quaternion_normalize(&integrator->q);
}

/**
 * 获取当前姿态四元数
 * 
 * @param integrator 积分器结构体指针
 * @param w, x, y, z 输出参数，四元数分量
 */
void angular_get_quaternion(const AngularIntegrator* integrator,
                           float* w, float* x, float* y, float* z) {
    *w = integrator->q.w;
    *x = integrator->q.x;
    *y = integrator->q.y;
    *z = integrator->q.z;
}

/**
 * 将四元数转换为欧拉角（Roll-Pitch-Yaw）
 * 
 * 欧拉角定义（ZYX顺序，即Yaw-Pitch-Roll）：
 * - Roll  (φ): 绕X轴旋转，范围 [-π, π]
 * - Pitch (θ): 绕Y轴旋转，范围 [-π/2, π/2]
 * - Yaw   (ψ): 绕Z轴旋转，范围 [-π, π]
 * 
 * @param q 四元数指针
 * @param roll, pitch, yaw 输出参数，欧拉角（单位：弧度）
 * 
 * 【万向锁问题】
 * 当pitch接近±90°时，会出现万向锁奇异性
 * 此时roll和yaw的定义不再独立
 */
void quaternion_to_euler(const Quaternion* q, 
                        float* roll, float* pitch, float* yaw) {
    // Roll (绕X轴)
    float sinr_cosp = 2.0f * (q->w * q->x + q->y * q->z);
    float cosr_cosp = 1.0f - 2.0f * (q->x * q->x + q->y * q->y);
    *roll = atan2f(sinr_cosp, cosr_cosp);
    
    // Pitch (绕Y轴)
    float sinp = 2.0f * (q->w * q->y - q->z * q->x);
    if (fabsf(sinp) >= 1.0f) {
        // 万向锁：pitch = ±90°
        *pitch = copysignf(M_PI / 2.0f, sinp);
    } else {
        *pitch = asinf(sinp);
    }
    
    // Yaw (绕Z轴)
    float siny_cosp = 2.0f * (q->w * q->z + q->x * q->y);
    float cosy_cosp = 1.0f - 2.0f * (q->y * q->y + q->z * q->z);
    *yaw = atan2f(siny_cosp, cosy_cosp);
}

/*===========================================================================
 * 使用示例
 *===========================================================================*/

/*

// 示例1：基本使用
void example_basic() {
    AngularIntegrator integrator;
    
    // 初始化：100Hz采样率，使用中点法
    angular_init(&integrator, 0.01f, 1);
    
    // 主循环
    while (1) {
        // 读取陀螺仪数据（rad/s）
        float gx, gy, gz;
        read_gyroscope(&gx, &gy, &gz);
        
        // 更新姿态
        angular_update(&integrator, gx, gy, gz);
        
        // 获取四元数
        float w, x, y, z;
        angular_get_quaternion(&integrator, &w, &x, &y, &z);
        
        // 转换为欧拉角（度）
        float roll, pitch, yaw;
        quaternion_to_euler(&integrator.q, &roll, &pitch, &yaw);
        roll *= 57.2958f;   // 弧度转角度
        pitch *= 57.2958f;
        yaw *= 57.2958f;
        
        delay(10);  // 10ms = 100Hz
    }
}

// 示例2：高精度应用
void example_high_accuracy() {
    AngularIntegrator integrator;
    
    // 使用RK4方法获得最高精度
    angular_init(&integrator, 0.001f, 2);  // 1000Hz, RK4
    
    int counter = 0;
    while (1) {
        float gx, gy, gz;
        read_gyroscope(&gx, &gy, &gz);
        
        angular_update(&integrator, gx, gy, gz);
        
        // 每1000次迭代额外归一化一次
        if (++counter >= 1000) {
            quaternion_normalize(&integrator.q);
            counter = 0;
        }
        
        delay(1);  // 1ms = 1000Hz
    }
}

// 示例3：与其他传感器融合
void example_sensor_fusion() {
    AngularIntegrator integrator;
    angular_init(&integrator, 0.01f, 1);
    
    while (1) {
        // 陀螺仪积分（快速跟踪）
        float gx, gy, gz;
        read_gyroscope(&gx, &gy, &gz);
        angular_update(&integrator, gx, gy, gz);
        
        // 每100次使用加速度计校正一次（避免漂移）
        static int acc_counter = 0;
        if (++acc_counter >= 100) {
            float ax, ay, az;
            read_accelerometer(&ax, &ay, &az);
            
            // 从加速度计估计roll和pitch
            float acc_roll = atan2f(ay, az);
            float acc_pitch = atan2f(-ax, sqrtf(ay*ay + az*az));
            
            // 简单互补滤波：98%陀螺仪 + 2%加速度计
            float roll, pitch, yaw;
            quaternion_to_euler(&integrator.q, &roll, &pitch, &yaw);
            
            roll = 0.98f * roll + 0.02f * acc_roll;
            pitch = 0.98f * pitch + 0.02f * acc_pitch;
            
            // 转回四元数并设置
            // （这里简化了，实际需要完整的欧拉角转四元数函数）
            
            acc_counter = 0;
        }
        
        delay(10);
    }
}

*/

/*******************************************************************************
 * Madgwick姿态解算算法 - C语言实现
 * 
 * 算法简介：
 * ==============================================================================
 * Madgwick算法是由Sebastian Madgwick提出的一种高效的姿态解算算法，适用于
 * IMU(惯性测量单元)和MARG(磁力、角速率、重力)传感器阵列。
 * 
 * 算法特点：
 * 1. 使用四元数表示姿态，避免了欧拉角的万向节死锁问题
 * 2. 采用梯度下降法优化姿态估计，计算效率高
 * 3. 只有一个可调参数(增益β)，调参简单
 * 4. 可在低采样率下运行，适合嵌入式系统
 * 5. 内置陀螺仪零漂补偿和磁干扰补偿
 * 
 * 算法原理：
 * ==============================================================================
 * 
 * 1. 基于角速度的姿态更新（预测步骤）
 * -----------------------------------------
 * 使用陀螺仪测量的角速度ω通过数值积分更新姿态：
 * 
 *   q_ω(t) = q(t-1) + q̇_ω(t)·Δt
 *          = q(t-1) + 0.5·q(t-1)⊗[0, ωx, ωy, ωz]·Δt
 * 
 * 其中：
 *   - q是四元数[qw, qx, qy, qz]
 *   - ω是三轴角速度[ωx, ωy, ωz] (单位: rad/s)
 *   - ⊗表示四元数乘法(Hamilton乘积)
 *   - Δt是采样周期
 * 
 * 2. 基于加速度计和磁力计的梯度下降优化（校正步骤）
 * ---------------------------------------------------------
 * 通过最小化目标函数来修正陀螺仪的累积误差：
 * 
 *   f(q, d, s) = q*·d·q - s
 * 
 * 其中：
 *   - d是地球坐标系的参考向量（重力或磁场）
 *   - s是传感器坐标系的测量向量
 *   - q*是q的共轭四元数
 * 
 * 使用梯度下降法求解：
 * 
 *   q_∇(t) = q(t-1) - β·(∇f / |∇f|)
 * 
 * 其中：
 *   - β是算法增益，控制收敛速度
 *   - ∇f = J^T·f 是目标函数的梯度
 *   - J是f对q的雅可比矩阵
 * 
 * 3. 融合估计
 * ------------
 * 最终姿态估计结合陀螺仪积分和梯度下降修正：
 * 
 *   q(t) = q(t-1) + [q̇_ω(t) - β·q̇_ε(t)]·Δt
 * 
 * 其中 q̇_ε = ∇f/|∇f| 是误差方向
 * 
 * 
 * IMU模式（仅使用陀螺仪和加速度计）
 * ==============================================================================
 * 
 * 目标函数（基于重力方向）：
 *   假设地球坐标系重力方向为 g = [0, 0, 1]
 *   加速度计归一化测量值为 a = [ax, ay, az]/|a|
 * 
 *   f_g(q, a) = [2(qx·qz - qw·qy) - ax]
 *               [2(qw·qx + qy·qz) - ay]
 *               [2(0.5 - qx² - qy²) - az]
 * 
 * 雅可比矩阵：
 *   J_g(q) = [-2qy    2qz   -2qw    2qx  ]
 *            [ 2qx    2qw    2qz    2qy  ]
 *            [ 0     -4qx   -4qy    0    ]
 * 
 * 梯度：∇f = J_g^T · f_g
 * 
 * 
 * MARG模式（使用陀螺仪、加速度计和磁力计）
 * ==============================================================================
 * 
 * 在IMU基础上增加磁场方向约束：
 *   地球磁场参考方向 b = [bx, 0, bz]（简化模型，忽略东向分量）
 *   磁力计归一化测量值 m = [mx, my, mz]/|m|
 * 
 * 目标函数扩展：
 *   f_{g,b}(q, a, b, m) = [f_g(q, a)]
 *                         [f_b(q, b, m)]
 * 
 * 其中 f_b 类似 f_g，但使用磁场向量
 * 
 * 
 * 参数调整指南：
 * ==============================================================================
 * 
 * 增益参数 β：
 *   - β ≈ √(3/4) · ω̄_β，其中 ω̄_β 是陀螺仪误差的平均值
 *   - IMU模式推荐值：0.033 (适合慢速运动)
 *   - MARG模式推荐值：0.041 (适合有磁干扰环境)
 *   - β越大：收敛越快，但对噪声更敏感
 *   - β越小：收敛越慢，但对噪声鲁棒性更好
 * 
 * 
 * 使用示例：
 * ==============================================================================
 * 
 * // 定义四元数和传感器数据
 * float q[4] = {1.0f, 0.0f, 0.0f, 0.0f};  // 初始姿态（单位四元数）
 * float gyr[3] = {wx, wy, wz};             // 陀螺仪数据 (rad/s)
 * float acc[3] = {ax, ay, az};             // 加速度计数据 (m/s²)
 * float mag[3] = {mx, my, mz};             // 磁力计数据 (可选)
 * float dt = 0.01f;                        // 采样周期 (100Hz)
 * float beta = 0.033f;                     // 算法增益
 * 
 * // IMU更新（仅陀螺仪+加速度计）
 * madgwick_update_imu(q, gyr, acc, beta, dt);
 * 
 * // MARG更新（陀螺仪+加速度计+磁力计）
 * madgwick_update_marg(q, gyr, acc, mag, beta, dt);
 * 
 * 
 * 注意事项：
 * ==============================================================================
 * 1. 输入的四元数必须是单位四元数（模为1）
 * 2. 加速度计数据会在函数内部自动归一化
 * 3. 磁力计数据会在函数内部自动归一化
 * 4. 建议在静止状态下初始化姿态
 * 5. 采样频率建议 ≥100Hz 以获得较好效果
 * 
 ******************************************************************************/

#include <math.h>

/*===========================================================================*/
/*                            辅助数学函数                                    */
/*===========================================================================*/

/**
 * @brief 向量归一化
 * @param v 输入输出向量，归一化后结果存回该数组
 * @param n 向量维度
 */
static void normalize(float *v, int n) {
    float norm = 0.0f;
    int i;
    
    // 计算向量的模
    for (i = 0; i < n; i++) {
        norm += v[i] * v[i];
    }
    norm = sqrtf(norm);
    
    // 避免除零
    if (norm > 0.0f) {
        for (i = 0; i < n; i++) {
            v[i] /= norm;
        }
    }
}

/**
 * @brief 四元数乘法（Hamilton乘积）
 * 
 * 四元数乘法定义：
 *   p⊗q = [pw·qw - px·qx - py·qy - pz·qz]
 *         [pw·qx + px·qw + py·qz - pz·qy]
 *         [pw·qy - px·qz + py·qw + pz·qx]
 *         [pw·qz + px·qy - py·qx + pz·qw]
 * 
 * @param p 第一个四元数 [pw, px, py, pz]
 * @param q 第二个四元数 [qw, qx, qy, qz]
 * @param result 结果四元数（输出）
 */
static void quaternion_multiply(const float *p, const float *q, float *result) {
    result[0] = p[0]*q[0] - p[1]*q[1] - p[2]*q[2] - p[3]*q[3];  // w分量
    result[1] = p[0]*q[1] + p[1]*q[0] + p[2]*q[3] - p[3]*q[2];  // x分量
    result[2] = p[0]*q[2] - p[1]*q[3] + p[2]*q[0] + p[3]*q[1];  // y分量
    result[3] = p[0]*q[3] + p[1]*q[2] - p[2]*q[1] + p[3]*q[0];  // z分量
}

/*===========================================================================*/
/*                      Madgwick IMU算法核心函数                              */
/*===========================================================================*/

/**
 * @brief Madgwick IMU姿态更新（仅使用陀螺仪和加速度计）
 * 
 * 算法流程：
 * 1. 从角速度计算四元数导数（预测）
 * 2. 使用加速度计测量计算梯度（校正）
 * 3. 融合预测和校正，更新四元数
 * 4. 归一化四元数
 * 
 * @param q 四元数 [qw, qx, qy, qz]（输入输出，会被更新）
 * @param gyr 陀螺仪测量 [ωx, ωy, ωz] 单位: rad/s
 * @param acc 加速度计测量 [ax, ay, az] 单位: m/s² (任意单位，会归一化)
 * @param beta 算法增益（推荐值 0.033）
 * @param dt 采样周期 (秒)
 */
void madgwick_update_imu(float *q, const float *gyr, const float *acc, float beta, float dt) {
    float qw = q[0], qx = q[1], qy = q[2], qz = q[3];
    float ax, ay, az;
    float gx = gyr[0], gy = gyr[1], gz = gyr[2];
    
    /* 步骤1: 从陀螺仪计算四元数变化率 q̇_ω */
    /* q̇ = 0.5 * q ⊗ [0, ωx, ωy, ωz] */
    float qDot[4];
    float omega_quat[4] = {0.0f, gx, gy, gz};  // 将角速度表示为纯四元数
    float temp_q[4] = {qw, qx, qy, qz};
    quaternion_multiply(temp_q, omega_quat, qDot);
    qDot[0] *= 0.5f;
    qDot[1] *= 0.5f;
    qDot[2] *= 0.5f;
    qDot[3] *= 0.5f;
    
    /* 步骤2: 使用加速度计校正（梯度下降） */
    /* 只有当加速度计有有效测量时才进行校正 */
    if ((acc[0] != 0.0f) || (acc[1] != 0.0f) || (acc[2] != 0.0f)) {
        float a_temp[3] = {acc[0], acc[1], acc[2]};
        normalize(a_temp, 3);  // 归一化加速度计测量
        ax = a_temp[0];
        ay = a_temp[1];
        az = a_temp[2];
        
        /* 计算目标函数 f_g(q, a) */
        /* 这是将重力向量从地球坐标系旋转到传感器坐标系，再与测量值对比 */
        float f[3];
        f[0] = 2.0f * (qx*qz - qw*qy) - ax;
        f[1] = 2.0f * (qw*qx + qy*qz) - ay;
        f[2] = 2.0f * (0.5f - qx*qx - qy*qy) - az;
        
        /* 计算雅可比矩阵 J_g 的转置与 f 的乘积（即梯度 ∇f） */
        /* J_g^T = [-2qy   2qx    0   ]
                   [ 2qz   2qw  -4qx  ]
                   [-2qw   2qz  -4qy  ]
                   [ 2qx   2qy   0   ] */
        float gradient[4];
        gradient[0] = -2.0f*qy*f[0] + 2.0f*qx*f[1];
        gradient[1] =  2.0f*qz*f[0] + 2.0f*qw*f[1] - 4.0f*qx*f[2];
        gradient[2] = -2.0f*qw*f[0] + 2.0f*qz*f[1] - 4.0f*qy*f[2];
        gradient[3] =  2.0f*qx*f[0] + 2.0f*qy*f[1];
        
        /* 归一化梯度 */
        normalize(gradient, 4);
        
        /* 应用梯度下降修正 */
        /* q̇ = q̇_ω - β·∇f */
        qDot[0] -= beta * gradient[0];
        qDot[1] -= beta * gradient[1];
        qDot[2] -= beta * gradient[2];
        qDot[3] -= beta * gradient[3];
    }
    
    /* 步骤3: 积分四元数变化率，更新姿态 */
    /* q(t) = q(t-1) + q̇·Δt */
    q[0] += qDot[0] * dt;
    q[1] += qDot[1] * dt;
    q[2] += qDot[2] * dt;
    q[3] += qDot[3] * dt;
    
    /* 步骤4: 归一化四元数（保证单位四元数约束） */
    normalize(q, 4);
}

/*===========================================================================*/
/*                     Madgwick MARG算法核心函数                             */
/*===========================================================================*/

/**
 * @brief Madgwick MARG姿态更新（使用陀螺仪、加速度计和磁力计）
 * 
 * 相比IMU版本，增加了磁场方向约束，可以估计完整的三维姿态（包括航向角）
 * 
 * 算法流程：
 * 1. 从角速度计算四元数导数（预测）
 * 2. 旋转磁力计测量到地球坐标系，确定参考磁场方向
 * 3. 使用加速度计和磁力计测量计算梯度（校正）
 * 4. 融合预测和校正，更新四元数
 * 5. 归一化四元数
 * 
 * @param q 四元数 [qw, qx, qy, qz]（输入输出，会被更新）
 * @param gyr 陀螺仪测量 [ωx, ωy, ωz] 单位: rad/s
 * @param acc 加速度计测量 [ax, ay, az] 单位: m/s²
 * @param mag 磁力计测量 [mx, my, mz] 单位: 任意（会归一化）
 * @param beta 算法增益（推荐值 0.041）
 * @param dt 采样周期 (秒)
 */
void madgwick_update_marg(float *q, const float *gyr, const float *acc, 
                          const float *mag, float beta, float dt) {
    float qw = q[0], qx = q[1], qy = q[2], qz = q[3];
    float ax, ay, az, mx, my, mz;
    float gx = gyr[0], gy = gyr[1], gz = gyr[2];
    
    /* 步骤1: 从陀螺仪计算四元数变化率 */
    float qDot[4];
    float omega_quat[4] = {0.0f, gx, gy, gz};
    float temp_q[4] = {qw, qx, qy, qz};
    quaternion_multiply(temp_q, omega_quat, qDot);
    qDot[0] *= 0.5f;
    qDot[1] *= 0.5f;
    qDot[2] *= 0.5f;
    qDot[3] *= 0.5f;
    
    /* 步骤2: 归一化传感器测量 */
    if ((acc[0] != 0.0f) || (acc[1] != 0.0f) || (acc[2] != 0.0f)) {
        float a_temp[3] = {acc[0], acc[1], acc[2]};
        float m_temp[3] = {mag[0], mag[1], mag[2]};
        normalize(a_temp, 3);
        normalize(m_temp, 3);
        ax = a_temp[0]; ay = a_temp[1]; az = a_temp[2];
        mx = m_temp[0]; my = m_temp[1]; mz = m_temp[2];
        
        /* 步骤3: 将磁力计测量旋转到地球坐标系 */
        /* h = q ⊗ [0, mx, my, mz] ⊗ q* */
        float h[4];
        float m_quat[4] = {0.0f, mx, my, mz};
        float q_conj[4] = {qw, -qx, -qy, -qz};  // 四元数共轭
        float temp1[4];
        quaternion_multiply(temp_q, m_quat, temp1);
        quaternion_multiply(temp1, q_conj, h);
        
        /* 步骤4: 计算参考磁场方向（只保留北向和下向分量，忽略东向） */
        /* b = [√(hx²+hy²), 0, hz] */
        float bx = sqrtf(h[1]*h[1] + h[2]*h[2]);
        float bz = h[3];
        
        /* 步骤5: 计算目标函数（包含重力和磁场约束） */
        float f[6];
        
        // 重力部分 f_g
        f[0] = 2.0f * (qx*qz - qw*qy) - ax;
        f[1] = 2.0f * (qw*qx + qy*qz) - ay;
        f[2] = 2.0f * (0.5f - qx*qx - qy*qy) - az;
        
        // 磁场部分 f_b
        f[3] = 2.0f*bx*(0.5f - qy*qy - qz*qz) + 2.0f*bz*(qx*qz - qw*qy) - mx;
        f[4] = 2.0f*bx*(qx*qy - qw*qz) + 2.0f*bz*(qw*qx + qy*qz) - my;
        f[5] = 2.0f*bx*(qw*qy + qx*qz) + 2.0f*bz*(0.5f - qx*qx - qy*qy) - mz;
        
        /* 步骤6: 计算雅可比矩阵转置与目标函数的乘积（梯度） */
        float gradient[4];
        gradient[0] = -2.0f*qy*f[0] + 2.0f*qx*f[1]
                     -2.0f*bz*qy*f[3] + (-2.0f*bx*qz + 2.0f*bz*qx)*f[4] 
                     + 2.0f*bx*qy*f[5];
        
        gradient[1] = 2.0f*qz*f[0] + 2.0f*qw*f[1] - 4.0f*qx*f[2]
                     + 2.0f*bz*qz*f[3] + (2.0f*bx*qy + 2.0f*bz*qw)*f[4]
                     + (2.0f*bx*qz - 4.0f*bz*qx)*f[5];
        
        gradient[2] = -2.0f*qw*f[0] + 2.0f*qz*f[1] - 4.0f*qy*f[2]
                     + (-4.0f*bx*qy - 2.0f*bz*qw)*f[3]
                     + (2.0f*bx*qx + 2.0f*bz*qz)*f[4]
                     + (2.0f*bx*qw - 4.0f*bz*qy)*f[5];
        
        gradient[3] = 2.0f*qx*f[0] + 2.0f*qy*f[1]
                     + (-4.0f*bx*qz + 2.0f*bz*qx)*f[3]
                     + (-2.0f*bx*qw + 2.0f*bz*qy)*f[4]
                     + 2.0f*bx*qx*f[5];
        
        /* 归一化梯度 */
        normalize(gradient, 4);
        
        /* 步骤7: 应用梯度修正 */
        qDot[0] -= beta * gradient[0];
        qDot[1] -= beta * gradient[1];
        qDot[2] -= beta * gradient[2];
        qDot[3] -= beta * gradient[3];
    }
    
    /* 步骤8: 积分并更新四元数 */
    q[0] += qDot[0] * dt;
    q[1] += qDot[1] * dt;
    q[2] += qDot[2] * dt;
    q[3] += qDot[3] * dt;
    
    /* 步骤9: 归一化四元数 */
    normalize(q, 4);
}

/*===========================================================================*/
/*                            姿态转换辅助函数                                */
/*===========================================================================*/

/**
 * @brief 将四元数转换为欧拉角（Roll-Pitch-Yaw）
 * 
 * 欧拉角定义（ZYX旋转顺序）：
 *   Roll (φ):  绕X轴旋转，范围 [-π, π]
 *   Pitch (θ): 绕Y轴旋转，范围 [-π/2, π/2]
 *   Yaw (ψ):   绕Z轴旋转，范围 [-π, π]
 * 
 * @param q 输入四元数 [qw, qx, qy, qz]
 * @param rpy 输出欧拉角 [roll, pitch, yaw] 单位：弧度
 */
void quaternion_to_euler(const float *q, float *rpy) {
    float qw = q[0], qx = q[1], qy = q[2], qz = q[3];
    
    // Roll (φ)
    float sinr_cosp = 2.0f * (qw*qx + qy*qz);
    float cosr_cosp = 1.0f - 2.0f * (qx*qx + qy*qy);
    rpy[0] = atan2f(sinr_cosp, cosr_cosp);
    
    // Pitch (θ)
    float sinp = 2.0f * (qw*qy - qz*qx);
    if (fabsf(sinp) >= 1.0f)
        rpy[1] = copysignf(M_PI / 2.0f, sinp);  // ±90度处理
    else
        rpy[1] = asinf(sinp);
    
    // Yaw (ψ)
    float siny_cosp = 2.0f * (qw*qz + qx*qy);
    float cosy_cosp = 1.0f - 2.0f * (qy*qy + qz*qz);
    rpy[2] = atan2f(siny_cosp, cosy_cosp);
}

/*===========================================================================*/
/*                                END OF FILE                                */
/*===========================================================================*/

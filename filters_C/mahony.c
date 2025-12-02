/*******************************************************************************
 * Mahony姿态解算算法 - C语言实现
 * 
 * 算法简介：
 * ==============================================================================
 * Mahony算法是由Robert Mahony等人提出的基于SO(3)的非线性互补滤波器，
 * 专门用于惯性和磁性传感器的姿态估计。
 * 
 * 算法特点：
 * 1. 使用显式互补滤波器结构，物理意义清晰
 * 2. 具有比例-积分(PI)控制器，可在线估计陀螺仪零偏
 * 3. 计算效率高，适合实时嵌入式应用
 * 4. 两个增益参数(Kp和Ki)，调参灵活
 * 5. 对传感器噪声具有良好的鲁棒性
 * 
 * 
 * 算法原理：
 * ==============================================================================
 * 
 * 传感器模型
 * -----------
 * 
 * 陀螺仪模型：
 *   Ω_y = Ω + b + μ
 * 
 *   其中：
 *   - Ω: 真实角速度
 *   - b: 零偏（缓变）
 *   - μ: 测量噪声
 * 
 * 加速度计模型（准静态条件）：
 *   a ≈ -R^T·g
 * 
 *   其中：
 *   - R: 旋转矩阵
 *   - g: 重力加速度 [0, 0, 9.8]^T
 * 
 * 磁力计模型：
 *   m = R^T·h
 * 
 *   其中：
 *   - h: 地磁场向量
 * 
 * 
 * 互补滤波器结构
 * ---------------
 * 
 * 互补滤波器的基本思想是融合两种估计：
 * 
 * 1. 高频估计（陀螺仪积分）：
 *    - 优点：快速响应，短期精度高
 *    - 缺点：存在累积漂移
 * 
 * 2. 低频估计（加速度计/磁力计）：
 *    - 优点：长期稳定，无漂移
 *    - 缺点：受噪声和加速度干扰影响
 * 
 * 
 * 显式互补滤波器（Explicit Complementary Filter）
 * -----------------------------------------------
 * 
 * 姿态更新方程（四元数形式）：
 * 
 *   q̇ = 0.5·q⊗p(Ω_y - b̂ + Kp·ω_mes)
 * 
 *   ḃ = -Ki·ω_mes
 * 
 * 其中：
 *   - q: 姿态四元数
 *   - b̂: 估计的陀螺仪零偏
 *   - ω_mes: 测量修正项（由加速度计和磁力计提供）
 *   - Kp: 比例增益（控制收敛速度）
 *   - Ki: 积分增益（控制零偏估计速度）
 * 
 * 
 * 测量修正项计算
 * ---------------
 * 
 * IMU模式（仅加速度计）：
 * 
 *   ω_mes = Ka·(a × v_a)
 * 
 *   其中：
 *   - a: 归一化加速度计测量
 *   - v_a = R^T·[0, 0, 1]: 预测的重力方向
 *   - Ka: 加速度计权重
 * 
 * MARG模式（加速度计+磁力计）：
 * 
 *   ω_mes = Ka·(a × v_a) + Km·(m × v_m)
 * 
 *   其中：
 *   - m: 归一化磁力计测量
 *   - v_m = R^T·h_ref: 预测的磁场方向
 *   - Km: 磁力计权重
 * 
 * 
 * 误差校正的物理意义
 * -------------------
 * 
 * 叉乘 (a × v_a) 的含义：
 *   - 如果测量值a和预测值v_a完全对齐，叉乘为零，无需修正
 *   - 如果存在误差，叉乘产生垂直于两个向量的修正轴
 *   - 修正轴的模正比于sin(θ)，其中θ是两向量夹角
 * 
 * PI控制器作用：
 *   - 比例项(Kp): 快速响应姿态误差
 *   - 积分项(Ki): 慢慢消除陀螺仪零偏
 * 
 * 
 * 参数调整指南：
 * ==============================================================================
 * 
 * Kp（比例增益）：
 *   - 典型值：0.5 ~ 2.0
 *   - 推荐值：1.0（平衡响应速度和稳定性）
 *   - 越大：对传感器修正响应越快，但对噪声敏感
 *   - 越小：更平滑，但对快速运动响应慢
 * 
 * Ki（积分增益）：
 *   - 典型值：0.0 ~ 0.5
 *   - 推荐值：0.3（IMU）, 0.0（MARG - 磁干扰环境）
 *   - 功能：补偿陀螺仪零偏
 *   - Ki=0: 不估计零偏（适合低成本陀螺仪）
 *   - Ki>0: 在线估计零偏（需要足够的激励）
 * 
 * 传感器权重：
 *   - 高动态加速：降低Ka，增大Kp
 *   - 磁干扰环境：降低Km或设为0
 *   - 静态环境：Ka=Km=1.0
 * 
 * 
 * Mahony vs Madgwick：
 * ==============================================================================
 * 
 * Mahony算法：
 *   + 有零偏估计能力（Ki>0时）
 *   + 两个增益，调参更灵活
 *   + 物理意义更直观（PI控制）
 *   - 需要调两个参数
 * 
 * Madgwick算法：
 *   + 只有一个增益参数
 *   + 梯度下降，理论基础扎实
 *   - 没有显式零偏估计
 *   - 在某些情况下收敛较慢
 * 
 * 实际应用建议：
 *   - 如果陀螺仪零偏明显：选Mahony（设Ki>0）
 *   - 如果想简单调参：选Madgwick
 *   - 如果需要最佳性能：两个都试，选表现更好的
 * 
 * 
 * 使用示例：
 * ==============================================================================
 * 
 * // 初始化
 * float q[4] = {1.0f, 0.0f, 0.0f, 0.0f};      // 姿态四元数
 * float gyro_bias[3] = {0.0f, 0.0f, 0.0f};    // 陀螺仪零偏估计
 * float gyr[3], acc[3], mag[3];                // 传感器数据
 * float dt = 0.01f;                            // 采样周期(100Hz)
 * float Kp = 1.0f, Ki = 0.3f;                  // PI增益
 * 
 * // 主循环
 * while (运行) {
 *     读取传感器(gyr, acc, mag);
 *     
 *     // IMU更新
 *     mahony_update_imu(q, gyr, acc, gyro_bias, Kp, Ki, dt);
 *     
 *     // 或 MARG更新
 *     mahony_update_marg(q, gyr, acc, mag, gyro_bias, Kp, Ki, dt);
 * }
 * 
 * 
 * 注意事项：
 * ==============================================================================
 * 1. 需要提供gyro_bias数组存储零偏估计（即使Ki=0）
 * 2. 如果不需要零偏估计，可以将Ki设为0
 * 3. 建议在静止状态初始化姿态
 * 4. 动态运动时，可以临时降低Ki避免错误的零偏估计
 * 5. 采样频率建议 ≥100Hz
 * 
 ******************************************************************************/

#include <math.h>

/*===========================================================================*/
/*                            辅助数学函数                                    */
/*===========================================================================*/

/**
 * @brief 向量归一化
 * @param v 输入输出向量
 * @param n 向量维度
 */
static void normalize(float *v, int n) {
    float norm = 0.0f;
    int i;
    for (i = 0; i < n; i++) {
        norm += v[i] * v[i];
    }
    norm = sqrtf(norm);
    if (norm > 0.0f) {
        for (i = 0; i < n; i++) {
            v[i] /= norm;
        }
    }
}

/**
 * @brief 向量叉乘 result = a × b
 * 
 * 叉乘的几何意义：
 *   - 方向：垂直于a和b构成的平面，遵循右手定则
 *   - 大小：|a×b| = |a|·|b|·sin(θ)，其中θ是a和b的夹角
 *   - 用途：在姿态估计中表示旋转误差的修正轴
 * 
 * @param a 第一个向量
 * @param b 第二个向量
 * @param result 叉乘结果（输出）
 */
static void cross_product(const float *a, const float *b, float *result) {
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
}

/**
 * @brief 四元数转旋转矩阵（转置形式）
 * 
 * 计算R^T，用于将地球坐标系的向量转换到传感器坐标系
 * 
 * @param q 四元数 [qw, qx, qy, qz]
 * @param R_transpose 输出的旋转矩阵转置（3x3，按行存储）
 */
static void quaternion_to_rotation_transpose(const float *q, float *R_transpose) {
    float qw = q[0], qx = q[1], qy = q[2], qz = q[3];
    float qw2 = qw*qw, qx2 = qx*qx, qy2 = qy*qy, qz2 = qz*qz;
    
    // R^T 第一行
    R_transpose[0] = qw2 + qx2 - qy2 - qz2;
    R_transpose[1] = 2.0f*(qx*qy + qw*qz);
    R_transpose[2] = 2.0f*(qx*qz - qw*qy);
    
    // R^T 第二行
    R_transpose[3] = 2.0f*(qx*qy - qw*qz);
    R_transpose[4] = qw2 - qx2 + qy2 - qz2;
    R_transpose[5] = 2.0f*(qy*qz + qw*qx);
    
    // R^T 第三行
    R_transpose[6] = 2.0f*(qx*qz + qw*qy);
    R_transpose[7] = 2.0f*(qy*qz - qw*qx);
    R_transpose[8] = qw2 - qx2 - qy2 + qz2;
}

/**
 * @brief 矩阵向量乘法 result = M * v
 * 
 * @param M 3x3矩阵（按行存储）
 * @param v 3维向量
 * @param result 输出向量
 */
static void matrix_vector_multiply(const float *M, const float *v, float *result) {
    result[0] = M[0]*v[0] + M[1]*v[1] + M[2]*v[2];
    result[1] = M[3]*v[0] + M[4]*v[1] + M[5]*v[2];
    result[2] = M[6]*v[0] + M[7]*v[1] + M[8]*v[2];
}

/**
 * @brief 四元数乘法 result = p ⊗ q
 * 
 * @param p 第一个四元数
 * @param q 第二个四元数
 * @param result 乘积四元数
 */
static void quaternion_multiply(const float *p, const float *q, float *result) {
    result[0] = p[0]*q[0] - p[1]*q[1] - p[2]*q[2] - p[3]*q[3];
    result[1] = p[0]*q[1] + p[1]*q[0] + p[2]*q[3] - p[3]*q[2];
    result[2] = p[0]*q[2] - p[1]*q[3] + p[2]*q[0] + p[3]*q[1];
    result[3] = p[0]*q[3] + p[1]*q[2] - p[2]*q[1] + p[3]*q[0];
}

/*===========================================================================*/
/*                       Mahony IMU算法核心函数                               */
/*===========================================================================*/

/**
 * @brief Mahony IMU姿态更新（仅使用陀螺仪和加速度计）
 * 
 * 算法流程：
 * 1. 计算预测的重力方向（基于当前姿态）
 * 2. 计算测量修正项（加速度计误差的叉乘）
 * 3. 使用PI控制器修正陀螺仪测量
 * 4. 更新零偏估计
 * 5. 积分四元数变化率
 * 6. 归一化四元数
 * 
 * @param q 四元数 [qw, qx, qy, qz]（输入输出）
 * @param gyr 陀螺仪测量 [ωx, ωy, ωz] (rad/s)
 * @param acc 加速度计测量 [ax, ay, az]
 * @param gyro_bias 陀螺仪零偏估计（输入输出）
 * @param Kp 比例增益（推荐值：1.0）
 * @param Ki 积分增益（推荐值：0.3，不需要零偏估计时设为0）
 * @param dt 采样周期 (秒)
 */
void mahony_update_imu(float *q, const float *gyr, const float *acc, 
                       float *gyro_bias, float Kp, float Ki, float dt) {
    float qw = q[0], qx = q[1], qy = q[2], qz = q[3];
    float ax, ay, az;
    float gx = gyr[0], gy = gyr[1], gz = gyr[2];
    
    /* 步骤1: 检查陀螺仪数据有效性 */
    float gyr_norm = sqrtf(gx*gx + gy*gy + gz*gz);
    if (gyr_norm == 0.0f) {
        return;  // 无角速度，不更新
    }
    
    /* 步骤2: 归一化加速度计测量 */
    float a_temp[3] = {acc[0], acc[1], acc[2]};
    float a_norm = sqrtf(a_temp[0]*a_temp[0] + a_temp[1]*a_temp[1] + a_temp[2]*a_temp[2]);
    
    if (a_norm > 0.0f) {
        /* 加速度计有效，进行姿态修正 */
        
        /* 归一化 */
        ax = a_temp[0] / a_norm;
        ay = a_temp[1] / a_norm;
        az = a_temp[2] / a_norm;
        
        /* 步骤3: 计算预测的重力方向（传感器坐标系） */
        /* v_a = R^T·[0, 0, 1] */
        /* 使用四元数计算（避免显式构造旋转矩阵）： */
        float v_a[3];
        v_a[0] = 2.0f * (qx*qz - qw*qy);        // R^T的第三列
        v_a[1] = 2.0f * (qy*qz + qw*qx);
        v_a[2] = qw*qw - qx*qx - qy*qy + qz*qz;
        
        /* 步骤4: 计算误差修正项 ω_mes = a × v_a */
        /* 这个叉乘给出了需要旋转的轴和大小 */
        float omega_mes[3];
        float a_vec[3] = {ax, ay, az};
        cross_product(a_vec, v_a, omega_mes);
        
        /* 步骤5: 更新零偏估计（积分项） */
        /* ḃ = -Ki·ω_mes */
        /* b(t) = b(t-1) + ḃ·Δt */
        if (Ki > 0.0f) {
            gyro_bias[0] -= Ki * omega_mes[0] * dt;
            gyro_bias[1] -= Ki * omega_mes[1] * dt;
            gyro_bias[2] -= Ki * omega_mes[2] * dt;
        }
        
        /* 步骤6: 应用PI修正到陀螺仪测量 */
        /* Ω = Ω_y - b̂ + Kp·ω_mes */
        gx = gx - gyro_bias[0] + Kp * omega_mes[0];
        gy = gy - gyro_bias[1] + Kp * omega_mes[1];
        gz = gz - gyro_bias[2] + Kp * omega_mes[2];
        
    } else {
        /* 加速度计无效，仅使用陀螺仪（去除零偏） */
        gx = gx - gyro_bias[0];
        gy = gy - gyro_bias[1];
        gz = gz - gyro_bias[2];
    }
    
    /* 步骤7: 计算四元数变化率 q̇ = 0.5·q⊗p(Ω) */
    float p[4] = {0.0f, gx, gy, gz};  // 纯四元数表示角速度
    float qDot[4];
    quaternion_multiply(q, p, qDot);
    qDot[0] *= 0.5f;
    qDot[1] *= 0.5f;
    qDot[2] *= 0.5f;
    qDot[3] *= 0.5f;
    
    /* 步骤8: 积分更新四元数 q(t) = q(t-1) + q̇·Δt */
    q[0] += qDot[0] * dt;
    q[1] += qDot[1] * dt;
    q[2] += qDot[2] * dt;
    q[3] += qDot[3] * dt;
    
    /* 步骤9: 归一化四元数 */
    normalize(q, 4);
}

/*===========================================================================*/
/*                      Mahony MARG算法核心函数                              */
/*===========================================================================*/

/**
 * @brief Mahony MARG姿态更新（使用陀螺仪、加速度计和磁力计）
 * 
 * 相比IMU版本，增加了磁场方向修正，可以估计完整的三维姿态
 * 
 * 算法流程：
 * 1. 计算预测的重力和磁场方向
 * 2. 计算测量修正项（加速度计和磁力计误差的叉乘和）
 * 3. 使用PI控制器修正陀螺仪测量
 * 4. 更新零偏估计
 * 5. 积分四元数变化率
 * 6. 归一化四元数
 * 
 * @param q 四元数 [qw, qx, qy, qz]（输入输出）
 * @param gyr 陀螺仪测量 [ωx, ωy, ωz] (rad/s)
 * @param acc 加速度计测量 [ax, ay, az]
 * @param mag 磁力计测量 [mx, my, mz]
 * @param gyro_bias 陀螺仪零偏估计（输入输出）
 * @param Kp 比例增益（推荐值：1.0）
 * @param Ki 积分增益（推荐值：0.0~0.3）
 * @param dt 采样周期 (秒)
 */
void mahony_update_marg(float *q, const float *gyr, const float *acc, 
                        const float *mag, float *gyro_bias, 
                        float Kp, float Ki, float dt) {
    float qw = q[0], qx = q[1], qy = q[2], qz = q[3];
    float ax, ay, az, mx, my, mz;
    float gx = gyr[0], gy = gyr[1], gz = gyr[2];
    
    /* 步骤1: 检查陀螺仪数据 */
    if ((gx == 0.0f) && (gy == 0.0f) && (gz == 0.0f)) {
        return;
    }
    
    /* 步骤2: 归一化加速度计和磁力计测量 */
    float a_temp[3] = {acc[0], acc[1], acc[2]};
    float m_temp[3] = {mag[0], mag[1], mag[2]};
    float a_norm = sqrtf(a_temp[0]*a_temp[0] + a_temp[1]*a_temp[1] + a_temp[2]*a_temp[2]);
    float m_norm = sqrtf(m_temp[0]*m_temp[0] + m_temp[1]*m_temp[1] + m_temp[2]*m_temp[2]);
    
    if (a_norm > 0.0f) {
        /* 加速度计有效 */
        ax = a_temp[0] / a_norm;
        ay = a_temp[1] / a_norm;
        az = a_temp[2] / a_norm;
        
        /* 如果磁力计无效，退化为IMU模式 */
        if (m_norm == 0.0f) {
            mahony_update_imu(q, gyr, acc, gyro_bias, Kp, Ki, dt);
            return;
        }
        
        /* 磁力计有效，继续MARG处理 */
        mx = m_temp[0] / m_norm;
        my = m_temp[1] / m_norm;
        mz = m_temp[2] / m_norm;
        
        /* 步骤3: 计算预测的重力方向 v_a */
        float v_a[3];
        v_a[0] = 2.0f * (qx*qz - qw*qy);
        v_a[1] = 2.0f * (qy*qz + qw*qx);
        v_a[2] = qw*qw - qx*qx - qy*qy + qz*qz;
        
        /* 步骤4: 将磁力计测量旋转到地球坐标系 */
        /* h = R·m */
        float R_transpose[9];
        quaternion_to_rotation_transpose(q, R_transpose);
        float h[3];
        matrix_vector_multiply(R_transpose, m_temp, h);  // h = R^T^T·m = R·m
        
        /* 步骤5: 计算参考磁场方向（归一化，只保留水平和垂直分量） */
        /* 地磁场参考：b = [√(hx²+hy²), 0, hz] */
        float bx = sqrtf(h[0]*h[0] + h[1]*h[1]);
        float bz = h[2];
        
        /* 步骤6: 计算预测的磁场方向 v_m = R^T·[bx, 0, bz] */
        float b_ref[3] = {bx, 0.0f, bz};
        float v_m[3];
        /* 直接计算R^T·b_ref */
        v_m[0] = R_transpose[0]*bx + R_transpose[2]*bz;
        v_m[1] = R_transpose[3]*bx + R_transpose[5]*bz;
        v_m[2] = R_transpose[6]*bx + R_transpose[8]*bz;
        
        /* 步骤7: 计算综合误差修正项 */
        /* ω_mes = (a × v_a) + (m × v_m) */
        float omega_mes_a[3], omega_mes_m[3], omega_mes[3];
        float a_vec[3] = {ax, ay, az};
        float m_vec[3] = {mx, my, mz};
        
        cross_product(a_vec, v_a, omega_mes_a);  // 加速度计修正
        cross_product(m_vec, v_m, omega_mes_m);  // 磁力计修正
        
        omega_mes[0] = omega_mes_a[0] + omega_mes_m[0];
        omega_mes[1] = omega_mes_a[1] + omega_mes_m[1];
        omega_mes[2] = omega_mes_a[2] + omega_mes_m[2];
        
        /* 步骤8: 更新零偏估计 */
        if (Ki > 0.0f) {
            gyro_bias[0] -= Ki * omega_mes[0] * dt;
            gyro_bias[1] -= Ki * omega_mes[1] * dt;
            gyro_bias[2] -= Ki * omega_mes[2] * dt;
        }
        
        /* 步骤9: 应用PI修正 */
        gx = gx - gyro_bias[0] + Kp * omega_mes[0];
        gy = gy - gyro_bias[1] + Kp * omega_mes[1];
        gz = gz - gyro_bias[2] + Kp * omega_mes[2];
        
    } else {
        /* 加速度计无效，仅去除零偏 */
        gx = gx - gyro_bias[0];
        gy = gy - gyro_bias[1];
        gz = gz - gyro_bias[2];
    }
    
    /* 步骤10: 计算四元数变化率并积分 */
    float p[4] = {0.0f, gx, gy, gz};
    float qDot[4];
    quaternion_multiply(q, p, qDot);
    qDot[0] *= 0.5f;
    qDot[1] *= 0.5f;
    qDot[2] *= 0.5f;
    qDot[3] *= 0.5f;
    
    q[0] += qDot[0] * dt;
    q[1] += qDot[1] * dt;
    q[2] += qDot[2] * dt;
    q[3] += qDot[3] * dt;
    
    /* 步骤11: 归一化 */
    normalize(q, 4);
}

/*===========================================================================*/
/*                                END OF FILE                                */
/*===========================================================================*/

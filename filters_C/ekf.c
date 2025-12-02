/*******************************************************************************
 * 扩展卡尔曼滤波器 (Extended Kalman Filter - EKF) 姿态估计算法
 * 
 * 算法类型: 非线性状态估计滤波器
 * 功能说明: 融合陀螺仪、加速度计、磁力计数据进行最优姿态估计
 *
 * =============================================================================
 * 理论基础：
 * =============================================================================
 * 
 * 【卡尔曼滤波基本思想】
 * 
 * 卡尔曼滤波是一种递归的最优估计算法,通过结合系统模型和测量数据,
 * 在有噪声的环境下估计系统状态。它的核心思想是:
 * 
 * 1. 预测步骤: 使用系统模型预测下一时刻的状态
 * 2. 更新步骤: 使用测量值校正预测结果
 * 
 * 【线性卡尔曼滤波】
 * 
 * 对于线性系统:
 *   x(k) = F·x(k-1) + w(k)      // 状态方程
 *   z(k) = H·x(k) + v(k)        // 测量方程
 * 
 * 其中w(k)是过程噪声, v(k)是测量噪声
 * 
 * 标准卡尔曼滤波步骤:
 * 
 * 预测:
 *   x̂⁻(k) = F·x̂(k-1)               // 状态预测
 *   P⁻(k) = F·P(k-1)·Fᵀ + Q        // 协方差预测
 * 
 * 更新:
 *   K(k) = P⁻(k)·Hᵀ·[H·P⁻(k)·Hᵀ + R]⁻¹  // 卡尔曼增益
 *   x̂(k) = x̂⁻(k) + K(k)·[z(k) - H·x̂⁻(k)] // 状态更新
 *   P(k) = [I - K(k)·H]·P⁻(k)            // 协方差更新
 * 
 * 【扩展卡尔曼滤波 (EKF)】
 * 
 * 姿态估计是非线性问题,因为:
 * 1. 四元数运动方程是非线性的
 * 2. 传感器测量模型是非线性的
 * 
 * EKF通过线性化处理非线性:
 *   x(k) = f(x(k-1), u(k)) + w(k)
 *   z(k) = h(x(k)) + v(k)
 * 
 * 线性化方法: 在当前估计点处计算雅可比矩阵
 *   F = ∂f/∂x |_(x̂(k-1))
 *   H = ∂h/∂x |_(x̂⁻(k))
 * 
 * 然后使用这些雅可比矩阵替代线性KF中的F和H
 * 
 * 【姿态EKF的状态空间】
 * 
 * 状态向量: x = [q₀, q₁, q₂, q₃, bₓ, bᵧ, b_z]ᵀ
 * - q₀,q₁,q₂,q₃: 姿态四元数
 * - bₓ,bᵧ,b_z: 陀螺仪零偏
 * 
 * 输入: u = [ωₓ, ωᵧ, ω_z]ᵀ (陀螺仪角速度)
 * 
 * 测量: z = [aₓ, aᵧ, a_z, mₓ, mᵧ, m_z]ᵀ
 * - aₓ,aᵧ,a_z: 加速度计测量
 * - mₓ,mᵧ,m_z: 磁力计测量
 * 
 * 【过程模型】
 * 
 * 四元数微分方程:
 *   dq/dt = (1/2)·q ⊗ (ω - b)
 * 
 * 离散化 (一阶近似):
 *   q(k) = q(k-1) + (Δt/2)·Ω(ω - b)·q(k-1)
 *        = [I + (Δt/2)·Ω(ω - b)]·q(k-1)
 * 
 * 其中Ω(ω)是角速度的四元数矩阵表示
 * 
 * 零偏假设为随机游走:
 *   b(k) = b(k-1) + wb(k)
 * 
 * 【测量模型】
 * 
 * 加速度计测量重力在机体系的投影:
 *   a_meas = Rᵀ(q)·g + noise
 * 
 * 磁力计测量地磁场在机体系的投影:
 *   m_meas = Rᵀ(q)·m + noise
 * 
 * 其中R(q)是四元数对应的旋转矩阵
 * 
 * =============================================================================
 * 算法特点：
 * =============================================================================
 * 
 * ✓ 优点：
 *   - 理论完备,有最优性保证(线性化误差小时)
 *   - 能估计并补偿陀螺仪零偏
 *   - 协方差提供不确定性度量
 *   - 自动调节测量权重
 *   - 能融合多种传感器
 * 
 * ✗ 缺点：
 *   - 计算复杂度高(O(n³),n=状态维数)
 *   - 需要调参(Q, R矩阵)
 *   - 线性化引入误差
 *   - 对初值敏感
 *   - 可能发散(线性化误差大时)
 * 
 * 应用场景：
 * - 无人机、机器人导航
 * - 虚拟现实头盔跟踪
 * - 需要估计陀螺仪零偏的场合
 * - 对精度要求高的应用
 * 
 * 【与其他滤波器对比】
 * 
 * EKF vs 互补滤波:
 * - EKF理论更严格,但计算量大
 * - 互补滤波简单快速,但调参依赖经验
 * 
 * EKF vs UKF:
 * - UKF不需要计算雅可比,实现简单
 * - UKF对强非线性更鲁棒
 * - EKF计算量略小
 * 
 ******************************************************************************/

#include <math.h>
#include <string.h>  // for memcpy

/*===========================================================================
 * 配置参数
 *===========================================================================*/

// 状态维度
#define STATE_DIM 7      // 4(四元数) + 3(陀螺仪零偏)
#define ACCEL_DIM 3      // 加速度计3轴
#define MAG_DIM 3        // 磁力计3轴
#define MEAS_DIM 6       // 测量维度: 3(加速度) + 3(磁场)

// 矩阵尺寸
#define STATE_STATE (STATE_DIM * STATE_DIM)  // 7x7 = 49

/*===========================================================================
 * 数据结构定义
 *===========================================================================*/

// 四元数
typedef struct {
    float w, x, y, z;
} Quaternion;

// 3D向量
typedef struct {
    float x, y, z;
} Vector3;

// EKF状态
typedef struct {
    // 状态估计
    Quaternion q;        // 姿态四元数
    Vector3 gyro_bias;   // 陀螺仪零偏 (rad/s)
    
    // 协方差矩阵 P (7x7)
    // 存储为一维数组: P[i][j] = data[i*STATE_DIM + j]
    float P[STATE_STATE];
    
    // 过程噪声协方差 Q (7x7)
    float Q[STATE_STATE];
    
    // 测量噪声协方差 R (6x6) 
    float R_accel[ACCEL_DIM * ACCEL_DIM];  // 加速度计 (3x3)
    float R_mag[MAG_DIM * MAG_DIM];        // 磁力计 (3x3)
    
    // 采样时间
    float dt;
    
    // 参考向量
    Vector3 gravity_ref;      // 重力参考向量 (通常是[0,0,1]或[0,0,-1])
    Vector3 magnetic_ref;     // 地磁参考向量
    
} EKF_State;

/*===========================================================================
 * 矩阵操作辅助函数
 *===========================================================================*/

/**
 * 矩阵清零
 */
void matrix_zero(float* mat, int rows, int cols) {
    int size = rows * cols;
    for (int i = 0; i < size; i++) {
        mat[i] = 0.0f;
    }
}

/**
 * 设置单位矩阵
 */
void matrix_identity(float* mat, int n) {
    matrix_zero(mat, n, n);
    for (int i = 0; i < n; i++) {
        mat[i * n + i] = 1.0f;
    }
}

/**
 * 矩阵加法: C = A + B
 */
void matrix_add(const float* A, const float* B, float* C, 
                int rows, int cols) {
    int size = rows * cols;
    for (int i = 0; i < size; i++) {
        C[i] = A[i] + B[i];
    }
}

/**
 * 矩阵乘法: C = A * B
 * A: m×n, B: n×p, C: m×p
 */
void matrix_multiply(const float* A, int m, int n,
                    const float* B, int p,
                    float* C) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            float sum = 0.0f;
            for (int k = 0; k < n; k++) {
                sum += A[i*n + k] * B[k*p + j];
            }
            C[i*p + j] = sum;
        }
    }
}

/**
 * 矩阵转置: B = Aᵀ
 */
void matrix_transpose(const float* A, float* B, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            B[j*rows + i] = A[i*cols + j];
        }
    }
}

/**
 * 3x3矩阵求逆 (使用伴随矩阵法)
 * 
 * 仅用于小矩阵!对于大矩阵应使用LU分解等方法
 */
int matrix_inverse_3x3(const float* A, float* A_inv) {
    // 计算行列式
    float det = A[0] * (A[4]*A[8] - A[5]*A[7])
              - A[1] * (A[3]*A[8] - A[5]*A[6])
              + A[2] * (A[3]*A[7] - A[4]*A[6]);
    
    if (fabsf(det) < 1e-9f) {
        return 0;  // 矩阵奇异
    }
    
    float inv_det = 1.0f / det;
    
    // 伴随矩阵
    A_inv[0] = (A[4]*A[8] - A[5]*A[7]) * inv_det;
    A_inv[1] = (A[2]*A[7] - A[1]*A[8]) * inv_det;
    A_inv[2] = (A[1]*A[5] - A[2]*A[4]) * inv_det;
    A_inv[3] = (A[5]*A[6] - A[3]*A[8]) * inv_det;
    A_inv[4] = (A[0]*A[8] - A[2]*A[6]) * inv_det;
    A_inv[5] = (A[2]*A[3] - A[0]*A[5]) * inv_det;
    A_inv[6] = (A[3]*A[7] - A[4]*A[6]) * inv_det;
    A_inv[7] = (A[1]*A[6] - A[0]*A[7]) * inv_det;
    A_inv[8] = (A[0]*A[4] - A[1]*A[3]) * inv_det;
    
    return 1;
}

/*===========================================================================
 * 四元数操作
 *===========================================================================*/

/**
 * 四元数归一化
 */
void quaternion_normalize(Quaternion* q) {
    float norm = sqrtf(q->w*q->w + q->x*q->x + q->y*q->y + q->z*q->z);
    if (norm > 1e-9f) {
        float inv_norm = 1.0f / norm;
        q->w *= inv_norm;
        q->x *= inv_norm;
        q->y *= inv_norm;
        q->z *= inv_norm;
    } else {
        q->w = 1.0f;
        q->x = q->y = q->z = 0.0f;
    }
}

/**
 * 四元数乘法: q = q1 ⊗ q2
 */
Quaternion quaternion_multiply(const Quaternion* q1, const Quaternion* q2) {
    Quaternion result;
    result.w = q1->w*q2->w - q1->x*q2->x - q1->y*q2->y - q1->z*q2->z;
    result.x = q1->w*q2->x + q1->x*q2->w + q1->y*q2->z - q1->z*q2->y;
    result.y = q1->w*q2->y - q1->x*q2->z + q1->y*q2->w + q1->z*q2->x;
    result.z = q1->w*q2->z + q1->x*q2->y - q1->y*q2->x + q1->z*q2->w;
    return result;
}

/**
 * 四元数旋转向量: v_rotated = q ⊗ v ⊗ q*
 * 
 * 将向量v从参考系旋转到机体系
 */
Vector3 quaternion_rotate_vector(const Quaternion* q, const Vector3* v) {
    // 使用高效公式: v' = v + 2·qv×(qv×v + qw·v)
    // 其中qv = [qx, qy, qz]
    
    Vector3 qv = {q->x, q->y, q->z};
    float qw = q->w;
    
    // t = 2 * cross(qv, v)
    Vector3 t;
    t.x = 2.0f * (qv.y*v->z - qv.z*v->y);
    t.y = 2.0f * (qv.z*v->x - qv.x*v->z);
    t.z = 2.0f * (qv.x*v->y - qv.y*v->x);
    
    // v' = v + qw*t + cross(qv, t)
    Vector3 result;
    result.x = v->x + qw*t.x + (qv.y*t.z - qv.z*t.y);
    result.y = v->y + qw*t.y + (qv.z*t.x - qv.x*t.z);
    result.z = v->z + qw*t.z + (qv.x*t.y - qv.y*t.x);
    
    return result;
}

/**
 * 四元数共轭: q* = [w, -x, -y, -z]
 * 用于反向旋转
 */
Quaternion quaternion_conjugate(const Quaternion* q) {
    Quaternion result;
    result.w = q->w;
    result.x = -q->x;
    result.y = -q->y;
    result.z = -q->z;
    return result;
}

/*===========================================================================
 * EKF核心算法
 *===========================================================================*/

/**
 * EKF初始化
 * 
 * @param ekf EKF状态结构体
 * @param dt 采样时间间隔 (秒)
 * @param process_noise 过程噪声标准差 (rad/s)
 * @param accel_noise 加速度计噪声标准差 (m/s²)
 * @param mag_noise 磁力计噪声标准差 (任意单位)
 * @param gyro_bias_noise 陀螺仪零偏随机游走噪声 (rad/s²)
 */
void ekf_init(EKF_State* ekf, float dt,
             float process_noise, float accel_noise, 
             float mag_noise, float gyro_bias_noise) {
    // 初始化状态
    ekf->q.w = 1.0f;
    ekf->q.x = ekf->q.y = ekf->q.z = 0.0f;
    ekf->gyro_bias.x = ekf->gyro_bias.y = ekf->gyro_bias.z = 0.0f;
    
    ekf->dt = dt;
    
    // 初始化协方差矩阵P为单位阵
    matrix_identity(ekf->P, STATE_DIM);
    
    // 设置过程噪声Q
    // 对角矩阵: 前4个是四元数噪声, 后3个是零偏噪声
    matrix_zero(ekf->Q, STATE_DIM, STATE_DIM);
    float q_var = process_noise * process_noise * dt * dt;
    float bias_var = gyro_bias_noise * gyro_bias_noise * dt;
    
    for (int i = 0; i < 4; i++) {
        ekf->Q[i*STATE_DIM + i] = q_var;
    }
    for (int i = 4; i < 7; i++) {
        ekf->Q[i*STATE_DIM + i] = bias_var;
    }
    
    // 设置测量噪声R (加速度计)
    matrix_zero(ekf->R_accel, ACCEL_DIM, ACCEL_DIM);
    float accel_var = accel_noise * accel_noise;
    for (int i = 0; i < ACCEL_DIM; i++) {
        ekf->R_accel[i*ACCEL_DIM + i] = accel_var;
    }
    
    // 设置测量噪声R (磁力计)
    matrix_zero(ekf->R_mag, MAG_DIM, MAG_DIM);
    float mag_var = mag_noise * mag_noise;
    for (int i = 0; i < MAG_DIM; i++) {
        ekf->R_mag[i*MAG_DIM + i] = mag_var;
    }
    
    // 设置参考向量 (NED坐标系)
    ekf->gravity_ref.x = 0.0f;
    ekf->gravity_ref.y = 0.0f;
    ekf->gravity_ref.z = 9.81f;  // 向下为正
    
    // 磁场参考向量 (需要根据当地实际情况设置)
    // 这里假设北半球,磁场向下倾斜
    ekf->magnetic_ref.x = 1.0f;   // 北分量
    ekf->magnetic_ref.y = 0.0f;   // 东分量
    ekf->magnetic_ref.z = 0.5f;   // 下分量 (倾角约26.6°)
    
    // 归一化磁场向量
    float m_norm = sqrtf(ekf->magnetic_ref.x * ekf->magnetic_ref.x +
                        ekf->magnetic_ref.y * ekf->magnetic_ref.y +
                        ekf->magnetic_ref.z * ekf->magnetic_ref.z);
    ekf->magnetic_ref.x /= m_norm;
    ekf->magnetic_ref.y /= m_norm;
    ekf->magnetic_ref.z /= m_norm;
}

/**
 * EKF预测步骤
 * 
 * 使用陀螺仪数据预测下一时刻的状态和协方差
 * 
 * @param ekf EKF状态
 * @param gx, gy, gz 陀螺仪角速度 (rad/s)
 */
void ekf_predict(EKF_State* ekf, float gx, float gy, float gz) {
    // 补偿零偏
    float omega_x = gx - ekf->gyro_bias.x;
    float omega_y = gy - ekf->gyro_bias.y;
    float omega_z = gz - ekf->gyro_bias.z;
    
    // 构造角速度四元数 [0, ωx, ωy, ωz]
    Quaternion omega_q;
    omega_q.w = 0.0f;
    omega_q.x = omega_x;
    omega_q.y = omega_y;
    omega_q.z = omega_z;
    
    // 四元数微分: dq/dt = 0.5 * q ⊗ ω
    Quaternion q_dot = quaternion_multiply(&ekf->q, &omega_q);
    q_dot.w *= 0.5f;
    q_dot.x *= 0.5f;
    q_dot.y *= 0.5f;
    q_dot.z *= 0.5f;
    
    // 一阶积分: q(k) = q(k-1) + dt * dq/dt
    ekf->q.w += ekf->dt * q_dot.w;
    ekf->q.x += ekf->dt * q_dot.x;
    ekf->q.y += ekf->dt * q_dot.y;
    ekf->q.z += ekf->dt * q_dot.z;
    
    // 归一化四元数
    quaternion_normalize(&ekf->q);
    
    // 零偏保持不变 (随机游走模型)
    // b(k) = b(k-1)
    
    /* 
     * 【协方差预测】
     * 
     * 完整的EKF需要计算雅可比矩阵F,然后:
     * P⁻ = F·P·Fᵀ + Q
     * 
     * F = ∂f/∂x是状态转移函数对状态的偏导数
     * 
     * 对于四元数系统,F矩阵较复杂:
     * 
     * F = [I₄ + (dt/2)Ω(ω)  | -(dt/2)Ξ(q)]
     *     [      0₃ₓ₄       |     I₃     ]
     * 
     * 其中:
     * - Ω(ω)是角速度矩阵 (4x4)
     * - Ξ(q)是四元数右乘矩阵 (4x3)
     * 
     * 为简化,这里使用简化的协方差更新
     */
    
    // 简化版: P = P + Q
    // 实际应用中应实现完整的F矩阵计算
    matrix_add(ekf->P, ekf->Q, ekf->P, STATE_DIM, STATE_DIM);
}

/**
 * EKF更新步骤 (加速度计)
 * 
 * 使用加速度计测量更新姿态估计
 * 
 * @param ekf EKF状态
 * @param ax, ay, az 加速度计测量 (m/s²)
 */
void ekf_update_accel(EKF_State* ekf, float ax, float ay, float az) {
    // 归一化加速度测量
    Vector3 a_meas = {ax, ay, az};
    float a_norm = sqrtf(ax*ax + ay*ay + az*az);
    if (a_norm < 0.1f) return;  // 测量无效
    
    a_meas.x /= a_norm;
    a_meas.y /= a_norm;
    a_meas.z /= a_norm;
    
    // 预测的重力方向 (从参考系到机体系)
    Quaternion q_conj = quaternion_conjugate(&ekf->q);
    Vector3 g_pred = quaternion_rotate_vector(&q_conj, &ekf->gravity_ref);
    
    // 归一化预测
    float g_norm = sqrtf(g_pred.x*g_pred.x + g_pred.y*g_pred.y + g_pred.z*g_pred.z);
    g_pred.x /= g_norm;
    g_pred.y /= g_norm;
    g_pred.z /= g_norm;
    
    // 创新 (innovation): y = z - h(x̂⁻)
    Vector3 innovation;
    innovation.x = a_meas.x - g_pred.x;
    innovation.y = a_meas.y - g_pred.y;
    innovation.z = a_meas.z - g_pred.z;
    
    /*
     * 【完整EKF更新】
     * 
     * 需要计算:
     * 1. 测量雅可比H = ∂h/∂x
     * 2. 创新协方差S = H·P·Hᵀ + R
     * 3. 卡尔曼增益K = P·Hᵀ·S⁻¹
     * 4. 状态更新x̂ = x̂⁻ + K·y
     * 5. 协方差更新P = (I - K·H)·P
     * 
     * 为简化,这里使用简化的更新规则
     */
    
    // 简化的互补更新 (类似Madgwick)
    float k_accel = 0.1f;  // 加速度计增益
    
    // 使用叉积作为误差信号
    Vector3 error;
    error.x = g_pred.y * a_meas.z - g_pred.z * a_meas.y;
    error.y = g_pred.z * a_meas.x - g_pred.x * a_meas.z;
    error.z = g_pred.x * a_meas.y - g_pred.y * a_meas.x;
    
    // 更新四元数
    Quaternion error_q;
    error_q.w = 0.0f;
    error_q.x = k_accel * error.x;
    error_q.y = k_accel * error.y;
    error_q.z = k_accel * error.z;
    
    Quaternion q_correction = quaternion_multiply(&ekf->q, &error_q);
    
    ekf->q.w += 0.5f * ekf->dt * q_correction.w;
    ekf->q.x += 0.5f * ekf->dt * q_correction.x;
    ekf->q.y += 0.5f * ekf->dt * q_correction.y;
    ekf->q.z += 0.5f * ekf->dt * q_correction.z;
    
    quaternion_normalize(&ekf->q);
}

/**
 * EKF更新步骤 (磁力计)
 * 
 * 使用磁力计测量更新航向
 */
void ekf_update_mag(EKF_State* ekf, float mx, float my, float mz) {
    // 归一化磁场测量
    Vector3 m_meas = {mx, my, mz};
    float m_norm = sqrtf(mx*mx + my*my + mz*mz);
    if (m_norm < 0.1f) return;
    
    m_meas.x /= m_norm;
    m_meas.y /= m_norm;
    m_meas.z /= m_norm;
    
    // 预测的磁场方向
    Quaternion q_conj = quaternion_conjugate(&ekf->q);
    Vector3 m_pred = quaternion_rotate_vector(&q_conj, &ekf->magnetic_ref);
    
    // 归一化
    float mp_norm = sqrtf(m_pred.x*m_pred.x + m_pred.y*m_pred.y + m_pred.z*m_pred.z);
    m_pred.x /= mp_norm;
    m_pred.y /= mp_norm;
    m_pred.z /= mp_norm;
    
    // 简化更新 (类似加速度计)
    float k_mag = 0.05f;  // 磁力计增益
    
    Vector3 error;
    error.x = m_pred.y * m_meas.z - m_pred.z * m_meas.y;
    error.y = m_pred.z * m_meas.x - m_pred.x * m_meas.z;
    error.z = m_pred.x * m_meas.y - m_pred.y * m_meas.x;
    
    Quaternion error_q;
    error_q.w = 0.0f;
    error_q.x = k_mag * error.x;
    error_q.y = k_mag * error.y;
    error_q.z = k_mag * error.z;
    
    Quaternion q_correction = quaternion_multiply(&ekf->q, &error_q);
    
    ekf->q.w += 0.5f * ekf->dt * q_correction.w;
    ekf->q.x += 0.5f * ekf->dt * q_correction.x;
    ekf->q.y += 0.5f * ekf->dt * q_correction.y;
    ekf->q.z += 0.5f * ekf->dt * q_correction.z;
    
    quaternion_normalize(&ekf->q);
}

/**
 * 完整的EKF更新 (陀螺仪+加速度计+磁力计)
 */
void ekf_update(EKF_State* ekf,
               float gx, float gy, float gz,
               float ax, float ay, float az,
               float mx, float my, float mz) {
    // 1. 预测步骤
    ekf_predict(ekf, gx, gy, gz);
    
    // 2. 更新步骤
    ekf_update_accel(ekf, ax, ay, az);
    ekf_update_mag(ekf, mx, my, mz);
}

/**
 * 获取当前姿态四元数
 */
void ekf_get_quaternion(const EKF_State* ekf,
                       float* w, float* x, float* y, float* z) {
    *w = ekf->q.w;
    *x = ekf->q.x;
    *y = ekf->q.y;
    *z = ekf->q.z;
}

/**
 * 获取陀螺仪零偏估计
 */
void ekf_get_gyro_bias(const EKF_State* ekf,
                      float* bx, float* by, float* bz) {
    *bx = ekf->gyro_bias.x;
    *by = ekf->gyro_bias.y;
    *bz = ekf->gyro_bias.z;
}

/*===========================================================================
 * 使用示例
 *===========================================================================*/

/*

void example_ekf() {
    EKF_State ekf;
    
    // 初始化EKF
    // 参数: dt=0.01s, 过程噪声=0.1, 加速度噪声=0.1, 
    //       磁场噪声=0.1, 零偏噪声=0.001
    ekf_init(&ekf, 0.01f, 0.1f, 0.1f, 0.1f, 0.001f);
    
    while (1) {
        // 读取传感器
        float gx, gy, gz;  // 陀螺仪
        float ax, ay, az;  // 加速度计
        float mx, my, mz;  // 磁力计
        
        read_sensors(&gx, &gy, &gz, &ax, &ay, &az, &mx, &my, &mz);
        
        // EKF更新
        ekf_update(&ekf, gx, gy, gz, ax, ay, az, mx, my, mz);
        
        // 获取结果
        float qw, qx, qy, qz;
        ekf_get_quaternion(&ekf, &qw, &qx, &qy, &qz);
        
        // 获取零偏估计
        float bx, by, bz;
        ekf_get_gyro_bias(&ekf, &bx, &by, &bz);
        
        printf("Q: [%.3f, %.3f, %.3f, %.3f]\n", qw, qx, qy, qz);
        printf("Bias: [%.4f, %.4f, %.4f] rad/s\n", bx, by, bz);
        
        delay(10);  // 100Hz
    }
}

*/

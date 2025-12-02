/*******************************************************************************
 * 倾角估计算法 (Tilt Estimation)
 * 
 * 算法类型: 静态姿态估计
 * 功能说明: 从加速度计和磁力计测量值直接计算姿态角
 *
 * =============================================================================
 * 理论基础：
 * =============================================================================
 * 
 * 倾角算法是最简单的姿态估计方法，适用于静止或低动态场景。
 * 它利用重力加速度和地磁场这两个已知的参考矢量来确定姿态。
 * 
 * 【基本假设】
 * 1. 加速度计只测量重力加速度（无外力加速度）
 * 2. 磁力计测量地磁场方向
 * 3. 传感器静止或低速运动
 * 
 * 【坐标系定义】
 * - 全局坐标系（NED）：North-East-Down
 *   - X轴：北
 *   - Y轴：东
 *   - Z轴：向下（重力方向）
 * 
 * - 机体坐标系：
 *   - X轴：前
 *   - Y轴：右
 *   - Z轴：下
 * 
 * 【Roll和Pitch计算】
 * 
 * 当传感器静止时，加速度计测量值 a = [ax, ay, az] 就是重力矢量。
 * 通过几何关系可以得到：
 * 
 * Roll  (φ) = atan2(ay, az)
 * Pitch (θ) = atan2(-ax, sqrt(ay² + az²))
 * 
 * 使用atan2而不是atan的优势：
 * - atan2考虑了所有四个象限
 * - 避免了除零问题
 * - 结果范围是[-π, π]
 * 
 * 【Yaw计算（需要磁力计）】
 * 
 * 由于重力只提供了两个自由度的信息（倾斜方向），
 * 无法确定绕重力轴的旋转（Yaw角）。
 * 
 * 需要使用磁力计来确定航向：
 * 1. 将磁力计读数旋转到水平面（补偿roll和pitch）
 * 2. 计算水平分量相对于磁北的角度
 * 
 * 倾斜补偿后的磁场分量：
 * bx = mx*cos(θ) + my*sin(φ)*sin(θ) + mz*cos(φ)*sin(θ)
 * by = my*cos(φ) - mz*sin(φ)
 * 
 * 航向角：
 * Yaw (ψ) = atan2(-by, bx)
 * 
 * =============================================================================
 * 算法特点：
 * =============================================================================
 * 
 * ✓ 优点：
 *   - 实现极简单，计算量极小
 *   - 不需要历史数据
 *   - 不会累积误差
 *   - 对传感器噪声鲁棒（可加低通滤波）
 * 
 * ✗ 缺点：
 *   - 只能用于静止或低动态场景
 *   - 运动加速度会产生大误差
 *   - 响应可能较慢（如果加滤波）
 *   - 磁场干扰影响航向
 * 
 * 应用场景：
 * - 水平仪、倾角仪
 * - 静态姿态监测
 * - 为其他滤波器提供初始姿态
 * - 低动态的航向测量
 * 
 * 【误差来源】
 * 1. 运动加速度：违反了"只有重力"的假设
 * 2. 磁场干扰：附近的铁磁物质、电流
 * 3. 传感器偏置：零偏、标定误差
 * 4. 传感器噪声：白噪声、温漂
 * 
 ******************************************************************************/

#include <math.h>

/*===========================================================================
 * 数据结构定义
 *===========================================================================*/

// 四元数结构
typedef struct {
    float w, x, y, z;
} Quaternion;

// 欧拉角结构（单位：弧度）
typedef struct {
    float roll;   // 横滚角 φ，绕X轴
    float pitch;  // 俯仰角 θ，绕Y轴
    float yaw;    // 偏航角 ψ，绕Z轴
} EulerAngles;

// Tilt估计器结构
typedef struct {
    EulerAngles angles;       // 当前欧拉角
    Quaternion quaternion;    // 当前四元数
    
    // 低通滤波参数（可选）
    float alpha;              // 滤波系数：0-1，越小越平滑
    int use_filter;           // 是否使用低通滤波
    
    // 磁偏角补偿（可选）
    float magnetic_declination;  // 磁偏角（弧度）
} TiltEstimator;

/*===========================================================================
 * 辅助数学函数
 *===========================================================================*/

/**
 * 向量归一化
 * 
 * 将向量缩放到单位长度
 * 归一化后的向量只保留方向信息，去除幅值影响
 */
void normalize_vector3(float* x, float* y, float* z) {
    float norm = sqrtf(*x * *x + *y * *y + *z * *z);
    if (norm > 1e-6f) {  // 避免除零
        float inv_norm = 1.0f / norm;
        *x *= inv_norm;
        *y *= inv_norm;
        *z *= inv_norm;
    }
}

/**
 * 低通滤波器（一阶IIR）
 * 
 * 公式：y[n] = α * x[n] + (1-α) * y[n-1]
 * 
 * α = 0: 完全平滑，无响应
 * α = 1: 无滤波，直接跟随
 * α = 0.1-0.3: 典型值，平衡响应和平滑
 * 
 * @param current 当前测量值
 * @param previous 上一次滤波输出
 * @param alpha 滤波系数
 * @return 滤波后的值
 */
float low_pass_filter(float current, float previous, float alpha) {
    return alpha * current + (1.0f - alpha) * previous;
}

/*===========================================================================
 * Roll和Pitch计算
 *===========================================================================*/

/**
 * 从加速度计计算Roll和Pitch角
 * 
 * 理论推导：
 * 
 * 在机体坐标系中，重力矢量g在NED坐标系中为[0, 0, g]
 * 经过旋转矩阵R后，在机体系中的投影为[ax, ay, az]
 * 
 * 根据旋转矩阵的几何意义：
 * - Roll角决定了重力在YZ平面的投影方向
 * - Pitch角决定了重力偏离Z轴的角度
 * 
 * 几何关系：
 *           ^
 *          / |
 *         /  | az
 *        /   |
 *       /θ___v
 *      ay
 * 
 * tan(φ) = ay / az  =>  φ = atan2(ay, az)
 * tan(θ) = -ax / sqrt(ay² + az²)  =>  θ = atan2(-ax, sqrt(ay² + az²))
 * 
 * 注意负号：
 * - 当传感器向前倾斜（机头向下），ax为正，但pitch应该为正
 * - 因此使用-ax
 * 
 * @param ax, ay, az 加速度计三轴读数（已归一化）
 * @param roll, pitch 输出的Roll和Pitch角（弧度）
 */
void calculate_roll_pitch(float ax, float ay, float az, 
                         float* roll, float* pitch) {
    // Roll角：绕X轴的旋转
    // 使用atan2保证结果在[-π, π]范围内
    *roll = atan2f(ay, az);
    
    // Pitch角：绕Y轴的旋转
    // 使用sqrt(ay²+az²)作为分母，表示重力在YZ平面的投影长度
    // 这样可以避免当roll=±90°时的奇异性
    float denominator = sqrtf(ay * ay + az * az);
    *pitch = atan2f(-ax, denominator);
    
    /* 
     * 【奇异性处理】
     * 
     * 当传感器完全倒置（az < 0且ay≈0）时，
     * atan2(ay, az)会在±π附近跳变
     * 
     * 如果需要更平滑的过渡，可以添加额外处理：
     * 
     * if (az < 0) {
     *     if (ay >= 0) {
     *         *roll = M_PI - *roll;
     *     } else {
     *         *roll = -M_PI - *roll;
     *     }
     * }
     */
}

/**
 * 另一种Roll/Pitch计算公式（避免三角函数）
 * 
 * 在某些低端MCU上，三角函数计算很慢
 * 可以使用近似公式：
 * 
 * roll ≈ ay / az  （当|roll|较小时）
 * pitch ≈ -ax     （当|pitch|较小时）
 * 
 * 这种方法只适用于小角度（<30°）的情况
 */
void calculate_roll_pitch_fast(float ax, float ay, float az,
                              float* roll, float* pitch) {
    // 快速近似（仅适用于小角度）
    if (fabsf(az) > 0.1f) {  // 避免除零
        *roll = ay / az;
        *pitch = -ax / az;
    } else {
        // 当az接近0时，使用精确计算
        *roll = atan2f(ay, az);
        *pitch = atan2f(-ax, sqrtf(ay*ay + az*az));
    }
}

/*===========================================================================
 * Yaw角计算（需要磁力计）
 *===========================================================================*/

/**
 * 计算倾斜补偿后的航向角（Yaw）
 * 
 * 步骤详解：
 * 
 * 1. 为什么需要倾斜补偿？
 *    当传感器倾斜时，磁力计测量的是地磁场在机体坐标系的投影
 *    如果直接用XY分量计算航向，会因为倾斜产生误差
 * 
 * 2. 倾斜补偿原理：
 *    通过roll和pitch角，将磁场矢量"旋转回"水平面
 *    相当于计算：如果传感器是水平放置的，磁场方向是什么
 * 
 * 3. 旋转矩阵推导：
 *    m_horizontal = Ry(-θ) * Rx(-φ) * m_body
 * 
 *    其中Rx和Ry是绕X轴和Y轴的旋转矩阵
 *    展开后得到补偿公式
 * 
 * 4. 航向计算：
 *    在水平面上，磁北方向就是X轴正方向
 *    航向 = atan2(-by, bx)
 *    负号是因为Y轴指向东（右手系）
 * 
 * @param mx, my, mz 磁力计三轴读数（已归一化）
 * @param roll, pitch 当前的Roll和Pitch角（弧度）
 * @return Yaw角（弧度）
 * 
 * 【磁偏角】
 * 地磁北极与地理北极不重合，产生磁偏角（Declination）
 * 磁偏角随地理位置变化：
 * - 东偏：磁北在地理北东侧（为正）
 * - 西偏：磁北在地理北西侧（为负）
 * 
 * 真北航向 = 磁北航向 + 磁偏角
 * 
 * 中国大陆磁偏角约为-5°到+10°
 */
float calculate_yaw(float mx, float my, float mz, 
                   float roll, float pitch) {
    // 预计算三角函数值（优化性能）
    float cos_roll = cosf(roll);
    float sin_roll = sinf(roll);
    float cos_pitch = cosf(pitch);
    float sin_pitch = sinf(pitch);
    
    /*
     * 倾斜补偿公式推导：
     * 
     * 绕X轴旋转-roll（补偿roll）：
     * [mx']   [1      0         0    ][mx]
     * [my'] = [0  cos(φ)   sin(φ)][my]
     * [mz']   [0 -sin(φ)   cos(φ)][mz]
     * 
     * 绕Y轴旋转-pitch（补偿pitch）：
     * [mx'']   [cos(θ)  0  -sin(θ)][mx']
     * [my''] = [   0     1     0    ][my']
     * [mz'']   [sin(θ)  0   cos(θ)][mz']
     * 
     * 组合后得到水平面的磁场分量：
     */
    
    // 水平面X分量（指向磁北）
    float bx = mx * cos_pitch + 
               my * sin_roll * sin_pitch + 
               mz * cos_roll * sin_pitch;
    
    // 水平面Y分量（指向磁东）
    float by = my * cos_roll - 
               mz * sin_roll;
    
    // 计算航向角
    // 注意：使用-by是因为我们要的是从北向东的角度
    float yaw = atan2f(-by, bx);
    
    return yaw;
    
    /*
     * 【航向角范围】
     * atan2返回值范围：[-π, π]
     * 
     * 对应关系：
     * - 0°:    北
     * - 90°:   东
     * - 180°:  南
     * - -90°:  西
     * 
     * 如果需要[0, 2π]范围：
     * if (yaw < 0) yaw += 2.0f * M_PI;
     * 
     * 如果需要[0°, 360°]：
     * yaw_deg = yaw * 180.0f / M_PI;
     * if (yaw_deg < 0) yaw_deg += 360.0f;
     */
}

/*===========================================================================
 * 欧拉角与四元数转换
 *===========================================================================*/

/**
 * 将Roll-Pitch-Yaw欧拉角转换为四元数
 * 
 * 转换公式（ZYX旋转顺序）：
 * 
 * 分别计算三个旋转的四元数，然后组合：
 * q = qz(ψ) ⊗ qy(θ) ⊗ qx(φ)
 * 
 * 使用半角公式简化计算
 * 
 * @param roll, pitch, yaw 欧拉角（弧度）
 * @param q 输出的四元数
 */
void euler_to_quaternion(float roll, float pitch, float yaw, Quaternion* q) {
    // 计算半角
    float cr = cosf(roll * 0.5f);
    float sr = sinf(roll * 0.5f);
    float cp = cosf(pitch * 0.5f);
    float sp = sinf(pitch * 0.5f);
    float cy = cosf(yaw * 0.5f);
    float sy = sinf(yaw * 0.5f);
    
    // 组合旋转
    q->w = cr * cp * cy + sr * sp * sy;
    q->x = sr * cp * cy - cr * sp * sy;
    q->y = cr * sp * cy + sr * cp * sy;
    q->z = cr * cp * sy - sr * sp * cy;
}

/*===========================================================================
 * 主要接口函数
 *===========================================================================*/

/**
 * 初始化Tilt估计器
 * 
 * @param estimator 估计器结构体指针
 * @param use_filter 是否使用低通滤波
 * @param filter_alpha 滤波系数（0-1），仅当use_filter=1时有效
 *                     推荐值：0.1-0.3
 * @param mag_declination 磁偏角（度），用于航向校正
 */
void tilt_init(TiltEstimator* estimator, 
              int use_filter, 
              float filter_alpha,
              float mag_declination) {
    // 初始化姿态为零
    estimator->angles.roll = 0.0f;
    estimator->angles.pitch = 0.0f;
    estimator->angles.yaw = 0.0f;
    
    estimator->quaternion.w = 1.0f;
    estimator->quaternion.x = 0.0f;
    estimator->quaternion.y = 0.0f;
    estimator->quaternion.z = 0.0f;
    
    // 滤波参数
    estimator->use_filter = use_filter;
    estimator->alpha = filter_alpha;
    
    // 磁偏角转换为弧度
    estimator->magnetic_declination = mag_declination * M_PI / 180.0f;
}

/**
 * 仅使用加速度计更新姿态（只有Roll和Pitch）
 * 
 * 适用场景：
 * - 不关心航向角
 * - 没有磁力计
 * - 磁场干扰严重
 * 
 * @param estimator 估计器结构体指针
 * @param ax, ay, az 加速度计三轴读数（任意单位，将被归一化）
 */
void tilt_update_imu(TiltEstimator* estimator, 
                    float ax, float ay, float az) {
    // 归一化加速度矢量
    normalize_vector3(&ax, &ay, &az);
    
    // 计算Roll和Pitch
    float roll, pitch;
    calculate_roll_pitch(ax, ay, az, &roll, &pitch);
    
    // 应用低通滤波（如果启用）
    if (estimator->use_filter) {
        roll = low_pass_filter(roll, estimator->angles.roll, 
                              estimator->alpha);
        pitch = low_pass_filter(pitch, estimator->angles.pitch, 
                               estimator->alpha);
    }
    
    // 更新角度
    estimator->angles.roll = roll;
    estimator->angles.pitch = pitch;
    // Yaw保持不变（或设为0）
    estimator->angles.yaw = 0.0f;
    
    // 转换为四元数
    euler_to_quaternion(roll, pitch, 0.0f, &estimator->quaternion);
}

/**
 * 使用加速度计和磁力计更新完整姿态（Roll-Pitch-Yaw）
 * 
 * 这是Tilt算法的完整版本
 * 
 * @param estimator 估计器结构体指针
 * @param ax, ay, az 加速度计三轴读数
 * @param mx, my, mz 磁力计三轴读数
 * 
 * 【使用注意】
 * 1. 加速度计和磁力计应该是在同一坐标系下的测量值
 * 2. 两个传感器的轴向应该对齐
 * 3. 磁力计需要校准（消除硬铁和软铁误差）
 */
void tilt_update_marg(TiltEstimator* estimator,
                     float ax, float ay, float az,
                     float mx, float my, float mz) {
    // 归一化加速度和磁场矢量
    normalize_vector3(&ax, &ay, &az);
    normalize_vector3(&mx, &my, &mz);
    
    // 计算Roll和Pitch
    float roll, pitch;
    calculate_roll_pitch(ax, ay, az, &roll, &pitch);
    
    // 计算Yaw（使用倾斜补偿）
    float yaw = calculate_yaw(mx, my, mz, roll, pitch);
    
    // 应用磁偏角补偿
    yaw += estimator->magnetic_declination;
    
    // 归一化yaw到[-π, π]
    while (yaw > M_PI) yaw -= 2.0f * M_PI;
    while (yaw < -M_PI) yaw += 2.0f * M_PI;
    
    // 应用低通滤波（如果启用）
    if (estimator->use_filter) {
        roll = low_pass_filter(roll, estimator->angles.roll, 
                              estimator->alpha);
        pitch = low_pass_filter(pitch, estimator->angles.pitch, 
                               estimator->alpha);
        
        // Yaw角需要特殊处理（避免±180°跳变）
        float yaw_diff = yaw - estimator->angles.yaw;
        // 归一化角度差到[-π, π]
        while (yaw_diff > M_PI) yaw_diff -= 2.0f * M_PI;
        while (yaw_diff < -M_PI) yaw_diff += 2.0f * M_PI;
        
        yaw = estimator->angles.yaw + 
              estimator->alpha * yaw_diff;
    }
    
    // 更新角度
    estimator->angles.roll = roll;
    estimator->angles.pitch = pitch;
    estimator->angles.yaw = yaw;
    
    // 转换为四元数
    euler_to_quaternion(roll, pitch, yaw, &estimator->quaternion);
}

/**
 * 获取欧拉角（角度）
 * 
 * @param estimator 估计器结构体指针
 * @param roll, pitch, yaw 输出参数，欧拉角（单位：度）
 */
void tilt_get_angles_deg(const TiltEstimator* estimator,
                        float* roll, float* pitch, float* yaw) {
    *roll = estimator->angles.roll * 180.0f / M_PI;
    *pitch = estimator->angles.pitch * 180.0f / M_PI;
    *yaw = estimator->angles.yaw * 180.0f / M_PI;
}

/**
 * 获取欧拉角（弧度）
 */
void tilt_get_angles_rad(const TiltEstimator* estimator,
                        float* roll, float* pitch, float* yaw) {
    *roll = estimator->angles.roll;
    *pitch = estimator->angles.pitch;
    *yaw = estimator->angles.yaw;
}

/**
 * 获取四元数
 */
void tilt_get_quaternion(const TiltEstimator* estimator,
                        float* w, float* x, float* y, float* z) {
    *w = estimator->quaternion.w;
    *x = estimator->quaternion.x;
    *y = estimator->quaternion.y;
    *z = estimator->quaternion.z;
}

/*===========================================================================
 * 使用示例
 *===========================================================================*/

/*

// 示例1：基本倾角测量（只有roll和pitch）
void example_basic_tilt() {
    TiltEstimator tilt;
    
    // 初始化：使用滤波，系数0.2，无磁偏角
    tilt_init(&tilt, 1, 0.2f, 0.0f);
    
    while (1) {
        // 读取加速度计
        float ax, ay, az;
        read_accelerometer(&ax, &ay, &az);
        
        // 更新姿态
        tilt_update_imu(&tilt, ax, ay, az);
        
        // 获取角度
        float roll, pitch, yaw;
        tilt_get_angles_deg(&tilt, &roll, &pitch, &yaw);
        
        printf("Roll: %.2f°, Pitch: %.2f°\n", roll, pitch);
        
        delay(50);  // 20Hz
    }
}

// 示例2：完整姿态（含航向）
void example_full_attitude() {
    TiltEstimator tilt;
    
    // 北京磁偏角约为-7°（西偏）
    tilt_init(&tilt, 1, 0.15f, -7.0f);
    
    while (1) {
        float ax, ay, az;
        float mx, my, mz;
        
        read_accelerometer(&ax, &ay, &az);
        read_magnetometer(&mx, &my, &mz);
        
        // 更新完整姿态
        tilt_update_marg(&tilt, ax, ay, az, mx, my, mz);
        
        // 获取角度
        float roll, pitch, yaw;
        tilt_get_angles_deg(&tilt, &roll, &pitch, &yaw);
        
        printf("Roll: %.2f°, Pitch: %.2f°, Yaw: %.2f°\n", 
               roll, pitch, yaw);
        
        delay(50);
    }
}

// 示例3：动态调整滤波
void example_adaptive_filter() {
    TiltEstimator tilt;
    tilt_init(&tilt, 1, 0.2f, 0.0f);
    
    while (1) {
        float ax, ay, az;
        read_accelerometer(&ax, &ay, &az);
        
        // 计算加速度幅值（检测运动）
        float acc_magnitude = sqrtf(ax*ax + ay*ay + az*az);
        
        // 接近1g时，传感器可能静止，增加滤波
        // 远离1g时，有额外加速度，减少滤波以保持响应
        if (fabsf(acc_magnitude - 9.81f) < 0.5f) {
            tilt.alpha = 0.1f;  // 强滤波
        } else {
            tilt.alpha = 0.5f;  // 弱滤波
        }
        
        tilt_update_imu(&tilt, ax, ay, az);
        
        float roll, pitch, yaw;
        tilt_get_angles_deg(&tilt, &roll, &pitch, &yaw);
        
        delay(20);
    }
}

*/

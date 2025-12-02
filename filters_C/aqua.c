/*******************************************************************************
 * AQUA算法 (Algebraic QUaternion Algorithm)
 * 
 * 算法类型: 代数解法姿态估计
 * 功能说明: 分别计算倾角四元数和航向四元数,然后组合
 * 
 * 理论基础:
 * AQUA将姿态分解为"倾斜"(tilt)和"航向"(heading)两部分独立计算,
 * 避免了磁干扰对roll/pitch的影响。
 * 
 * 优点: 抗磁干扰,计算速度快
 * 缺点: 需要两次四元数计算
 * 适用: 磁场干扰环境,无人机导航
 ******************************************************************************/

#include <math.h>

void aqua_update(const float* acc, const float* mag, 
                float mag_dip, float* q_out) {
    // 归一化加速度
    float ax=acc[0], ay=acc[1], az=acc[2];
    float a_norm = sqrtf(ax*ax + ay*ay + az*az);
    if(a_norm < 0.01f) return;
    ax /= a_norm; ay /= a_norm; az /= a_norm;
    
    // 计算倾斜四元数(从加速度)
    float roll = atan2f(ay, az);
    float pitch = atan2f(-ax, sqrtf(ay*ay+az*az));
    
    float cr = cosf(roll*0.5f), sr = sinf(roll*0.5f);
    float cp = cosf(pitch*0.5f), sp = sinf(pitch*0.5f);
    
    // 倾斜四元数
    float q_tilt[4] = {
        cr*cp,        // w
        sr*cp,        // x
        cr*sp,        // y  
        -sr*sp        // z
    };
    
    // 归一化磁场
    float mx=mag[0], my=mag[1], mz=mag[2];
    float m_norm = sqrtf(mx*mx + my*my + mz*mz);
    if(m_norm < 0.01f) {
        q_out[0]=q_tilt[0]; q_out[1]=q_tilt[1];
        q_out[2]=q_tilt[2]; q_out[3]=q_tilt[3];
        return;
    }
    mx /= m_norm; my /= m_norm; mz /= m_norm;
    
    // 倾斜补偿
    float cos_r=cosf(roll), sin_r=sinf(roll);
    float cos_p=cosf(pitch), sin_p=sinf(pitch);
    
    float bx = mx*cos_p + my*sin_r*sin_p + mz*cos_r*sin_p;
    float by = my*cos_r - mz*sin_r;
    
    // 计算航向
    float yaw = atan2f(-by, bx);
    float cy = cosf(yaw*0.5f), sy = sinf(yaw*0.5f);
    
    // 组合四元数: q = q_tilt * q_yaw
    q_out[0] = q_tilt[0]*cy - q_tilt[3]*sy;  // w
    q_out[1] = q_tilt[1]*cy + q_tilt[2]*sy;  // x
    q_out[2] = q_tilt[2]*cy - q_tilt[1]*sy;  // y
    q_out[3] = q_tilt[3]*cy + q_tilt[0]*sy;  // z
}

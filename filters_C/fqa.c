/*******************************************************************************
 * FQA算法 (Factored Quaternion Algorithm)
 * 
 * 算法类型: 分解四元数算法
 * 功能说明: 将姿态四元数分解为elevation、roll和azimuth三部分
 * 
 * 理论基础:
 * FQA与TRIAD等效,但输出四元数而非旋转矩阵。
 * 通过将姿态分解为三个独立旋转,避免了磁场干扰对倾角的影响。
 * 
 * 优点: 磁干扰只影响航向
 * 缺点: 不适用于高动态环境
 * 适用: 静态/准静态姿态测量
 ******************************************************************************/

#include <math.h>

void fqa_update(const float* acc, const float* mag, float* q_out) {
    // 归一化加速度
    float ax=acc[0], ay=acc[1], az=acc[2];
    float a_norm = sqrtf(ax*ax + ay*ay + az*az);
    if(a_norm < 0.01f) return;
    ax/=a_norm; ay/=a_norm; az/=a_norm;
    
    // 计算elevation(俯仰)四元数
    float sin_theta = ax;
    float cos_theta = sqrtf(1.0f - sin_theta*sin_theta);
    
    float sin_theta_2 = (sin_theta >= 0 ? 1.0f : -1.0f) * 
                        sqrtf((1.0f - cos_theta) * 0.5f);
    float cos_theta_2 = sqrtf((1.0f + cos_theta) * 0.5f);
    
    float q_e[4] = {cos_theta_2, 0.0f, sin_theta_2, 0.0f};
    
    // 计算roll四元数
    float sin_phi = 0.0f, cos_phi = 1.0f;
    if(fabsf(cos_theta) > 0.01f) {
        sin_phi = -ay / cos_theta;
        cos_phi = -az / cos_theta;
    }
    
    float phi = atan2f(sin_phi, cos_phi);
    float sin_phi_2 = sinf(phi * 0.5f);
    float cos_phi_2 = cosf(phi * 0.5f);
    
    float q_r[4] = {cos_phi_2, sin_phi_2, 0.0f, 0.0f};
    
    // 组合elevation和roll
    float q_er[4];
    q_er[0] = q_r[0]*q_e[0] - q_r[1]*q_e[1] - q_r[2]*q_e[2] - q_r[3]*q_e[3];
    q_er[1] = q_r[0]*q_e[1] + q_r[1]*q_e[0] + q_r[2]*q_e[3] - q_r[3]*q_e[2];
    q_er[2] = q_r[0]*q_e[2] - q_r[1]*q_e[3] + q_r[2]*q_e[0] + q_r[3]*q_e[1];
    q_er[3] = q_r[0]*q_e[3] + q_r[1]*q_e[2] - q_r[2]*q_e[1] + q_r[3]*q_e[0];
    
    // 如果有磁力计,计算azimuth
    if(mag != NULL) {
        float mx=mag[0], my=mag[1], mz=mag[2];
        float m_norm = sqrtf(mx*mx + my*my + mz*mz);
        if(m_norm > 0.01f) {
            mx/=m_norm; my/=m_norm; mz/=m_norm;
            
            // 简化的航向计算
            float yaw = atan2f(-my, mx);
            float cy = cosf(yaw*0.5f), sy = sinf(yaw*0.5f);
            
            // 组合最终四元数
            q_out[0] = q_er[0]*cy - q_er[3]*sy;
            q_out[1] = q_er[1]*cy + q_er[2]*sy;
            q_out[2] = q_er[2]*cy - q_er[1]*sy;
            q_out[3] = q_er[3]*cy + q_er[0]*sy;
            return;
        }
    }
    
    // 无磁力计,只输出倾角
    q_out[0]=q_er[0]; q_out[1]=q_er[1];
    q_out[2]=q_er[2]; q_out[3]=q_er[3];
}

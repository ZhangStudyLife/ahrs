/*******************************************************************************
 * SAAM算法 (Super-fast Attitude from Accelerometer and Magnetometer)
 * 
 * 算法类型: 极简化Wahba问题解法
 * 功能说明: 通过简化计算快速求解姿态四元数
 * 
 * 理论基础:
 * SAAM是对Davenport's q-method的极度简化,将特征值问题简化为
 * 几个浮点运算,适合计算资源极其受限的平台。
 * 
 * 优点: 计算量最小,速度极快
 * 缺点: 精度略低于完整方法
 * 适用: 8位单片机,实时性要求极高的场合
 ******************************************************************************/

#include <math.h>

void saam_update(const float* acc, const float* mag, float* q_out) {
    // 归一化输入
    float ax=acc[0], ay=acc[1], az=acc[2];
    float a_norm = sqrtf(ax*ax + ay*ay + az*az);
    if(a_norm < 0.01f) return;
    ax/=a_norm; ay/=a_norm; az/=a_norm;
    
    float mx=mag[0], my=mag[1], mz=mag[2];
    float m_norm = sqrtf(mx*mx + my*my + mz*mz);
    if(m_norm < 0.01f) return;
    mx/=m_norm; my/=m_norm; mz/=m_norm;
    
    // 计算参考磁场水平和垂直分量
    float alpha = ax*mx + ay*my + az*mz;  // 点积
    float m_N = sqrtf(1.0f - alpha*alpha);  // 水平分量
    float m_D = alpha;                       // 垂直分量
    
    // SAAM基本解(直接公式)
    q_out[0] = -ay*(m_N + mx) + ax*my;              // w
    q_out[1] = (az - 1.0f)*(m_N + mx) + ax*(m_D - mz);  // x
    q_out[2] = (az - 1.0f)*my + ay*(m_D - mz);          // y
    q_out[3] = az*m_D - ax*m_N - mz;                    // z
    
    // 归一化
    float norm = sqrtf(q_out[0]*q_out[0] + q_out[1]*q_out[1] +
                       q_out[2]*q_out[2] + q_out[3]*q_out[3]);
    if(norm > 1e-6f) {
        float inv = 1.0f / norm;
        q_out[0] *= inv; q_out[1] *= inv;
        q_out[2] *= inv; q_out[3] *= inv;
    }
}

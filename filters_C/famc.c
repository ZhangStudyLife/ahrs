/*******************************************************************************
 * FAMC算法 (Fast Accelerometer-Magnetometer Combination)
 * 
 * 算法类型: 快速ACC-MAG组合
 * 功能说明: 简化Wahba问题的解析解
 * 
 * 理论基础:
 * FAMC给出了ACC-MAG组合的解析特征值,避免了数值特征值分解,
 * 是Davenport's q-method的快速实现。
 * 
 * 优点: 速度快,精度高
 * 缺点: 仅限ACC-MAG组合
 * 适用: 静态姿态,快速计算场合
 ******************************************************************************/

#include <math.h>

void famc_update(const float* acc, const float* mag, float* q_out) {
    // 归一化
    float ax=acc[0], ay=acc[1], az=acc[2];
    float a_norm = sqrtf(ax*ax+ay*ay+az*az);
    if(a_norm<0.01f) return;
    ax/=a_norm; ay/=a_norm; az/=a_norm;
    
    float mx=mag[0], my=mag[1], mz=mag[2];
    float m_norm = sqrtf(mx*mx+my*my+mz*mz);
    if(m_norm<0.01f) return;
    mx/=m_norm; my/=m_norm; mz/=m_norm;
    
    // 参考向量(NED)
    float cos_dip = 0.866f;  // 30度磁倾角
    float sin_dip = 0.5f;
    
    // 构建B矩阵(简化)
    float B11 = 0.5f * (ax*0.0f + mx*cos_dip);
    float B13 = 0.5f * (ax*1.0f + mx*sin_dip);
    float B21 = 0.5f * (ay*0.0f + my*cos_dip);
    float B23 = 0.5f * (ay*1.0f + my*sin_dip);
    float B31 = 0.5f * (az*0.0f + mz*cos_dip);
    float B33 = 0.5f * (az*1.0f + mz*sin_dip);
    
    // 最大特征值为1(ACC-MAG组合的特殊性质)
    // 计算特征向量
    float a = B23;
    float b = B21;
    float c = B13;
    
    // 四元数(特征向量)
    q_out[0] = a;  // w
    q_out[1] = b;  // x
    q_out[2] = c;  // y
    q_out[3] = -1.0f;  // z
    
    // 归一化
    float norm = sqrtf(a*a + b*b + c*c + 1.0f);
    if(norm > 1e-6f) {
        float inv = 1.0f/norm;
        q_out[0]*=inv; q_out[1]*=inv;
        q_out[2]*=inv; q_out[3]*=inv;
    }
}

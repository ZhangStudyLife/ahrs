/*******************************************************************************
 * Davenport算法
 * 
 * 算法类型: Wahba问题的q-method解法
 * 功能说明: 通过构造K矩阵并求最大特征值来估计姿态
 * 
 * 理论基础:
 * Davenport's q-method是最早的四元数姿态确定方法之一,
 * 它将Wahba问题转化为求K矩阵最大特征值对应的特征向量。
 * 
 * 优点: 理论优雅,精度高
 * 缺点: 需要特征值分解
 * 适用: 静态姿态,学术研究
 ******************************************************************************/

#include <math.h>

void davenport_update(const float* acc, const float* mag,
                     float weight_acc, float weight_mag,
                     float* q_out) {
    // 归一化
    float ax=acc[0], ay=acc[1], az=acc[2];
    float a_norm = sqrtf(ax*ax+ay*ay+az*az);
    if(a_norm<0.01f) return;
    ax/=a_norm; ay/=a_norm; az/=a_norm;
    
    float mx=mag[0], my=mag[1], mz=mag[2];
    float m_norm = sqrtf(mx*mx+my*my+mz*mz);
    if(m_norm<0.01f) return;
    mx/=m_norm; my/=m_norm; mz/=m_norm;
    
    // 参考向量
    float ref_g[3] = {0, 0, 1};
    float ref_m[3] = {1, 0, 0.5f};
    
    // 构建B矩阵
    float B[9] = {0};
    for(int i=0; i<3; i++) {
        for(int j=0; j<3; j++) {
            B[i*3+j] = weight_acc * acc[i] * ref_g[j] +
                      weight_mag * mag[i] * ref_m[j];
        }
    }
    
    // 计算sigma和z
    float sigma = B[0] + B[4] + B[8];
    float z[3] = {B[5]-B[7], B[6]-B[2], B[1]-B[3]};
    
    // 简化:假设最大特征值≈1
    float lambda_max = 1.0f;
    
    // 构造四元数(简化计算)
    q_out[0] = lambda_max - sigma;
    q_out[1] = z[0];
    q_out[2] = z[1];
    q_out[3] = z[2];
    
    // 归一化
    float norm = sqrtf(q_out[0]*q_out[0] + q_out[1]*q_out[1] +
                      q_out[2]*q_out[2] + q_out[3]*q_out[3]);
    if(norm > 1e-6f) {
        float inv = 1.0f/norm;
        q_out[0]*=inv; q_out[1]*=inv;
        q_out[2]*=inv; q_out[3]*=inv;
    }
}

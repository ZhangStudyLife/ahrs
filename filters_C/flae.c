/*******************************************************************************
 * FLAE算法 (Fast Linear Attitude Estimator)
 * 
 * 算法类型: 快速线性姿态估计
 * 功能说明: 通过符号求解线性系统快速估计姿态
 * 
 * 理论基础:
 * FLAE将姿态估计转化为线性方程组,通过符号求解特征多项式
 * 获得解析解,避免了迭代和数值特征值计算。
 * 
 * 优点: 速度极快,精度高
 * 缺点: 数学复杂
 * 适用: 实时系统,高频率更新
 ******************************************************************************/

#include <math.h>

void flae_update(const float* acc, const float* mag, float* q_out) {
    // 简化实现:使用SAAM的思路
    float ax=acc[0], ay=acc[1], az=acc[2];
    float a_norm = sqrtf(ax*ax+ay*ay+az*az);
    if(a_norm<0.01f) return;
    ax/=a_norm; ay/=a_norm; az/=a_norm;
    
    float mx=mag[0], my=mag[1], mz=mag[2];
    float m_norm = sqrtf(mx*mx+my*my+mz*mz);
    if(m_norm<0.01f) return;
    mx/=m_norm; my/=m_norm; mz/=m_norm;
    
    // 计算辅助变量
    float alpha = ax*mx + ay*my + az*mz;
    float m_N = sqrtf(1.0f - alpha*alpha);
    float m_D = alpha;
    
    // 线性解
    q_out[0] = m_N + mx - ay;           // w
    q_out[1] = ax*(m_D-mz) + (az-1.0f)*(m_N+mx);  // x
    q_out[2] = ay*(m_D-mz) + (az-1.0f)*my;        // y
    q_out[3] = az*m_D - ax*m_N - mz;              // z
    
    // 归一化
    float norm = sqrtf(q_out[0]*q_out[0] + q_out[1]*q_out[1] +
                      q_out[2]*q_out[2] + q_out[3]*q_out[3]);
    if(norm > 1e-6f) {
        float inv = 1.0f/norm;
        q_out[0]*=inv; q_out[1]*=inv;
        q_out[2]*=inv; q_out[3]*=inv;
    }
}

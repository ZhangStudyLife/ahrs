/*******************************************************************************
 * OLEQ算法 (Optimal Linear Estimator of Quaternion)
 * 
 * 算法类型: 优化线性估计
 * 功能说明: 通过迭代矩阵旋转快速收敛到最优四元数
 * 
 * 理论基础:
 * OLEQ基于Brouwer不动点定理,通过不断应用旋转矩阵使随机初始
 * 四元数收敛到最优解。相比QUEST不需要特征值计算。
 * 
 * 优点: 无需特征值,实现简单
 * 缺点: 需要迭代(通常5-10次)
 * 适用: 静态姿态,QUEST的替代方案
 ******************************************************************************/

#include <math.h>

void oleq_update(const float* acc, const float* mag,
                float weight_acc, float weight_mag,
                float* q_out) {
    // 归一化输入
    float ax=acc[0], ay=acc[1], az=acc[2];
    float a_norm = sqrtf(ax*ax+ay*ay+az*az);
    if(a_norm<0.01f) return;
    ax/=a_norm; ay/=a_norm; az/=a_norm;
    
    float mx=mag[0], my=mag[1], mz=mag[2];
    float m_norm = sqrtf(mx*mx+my*my+mz*mz);
    if(m_norm<0.01f) return;
    mx/=m_norm; my/=m_norm; mz/=m_norm;
    
    // 参考向量(NED)
    float ref_g[3] = {0, 0, 1};    // 重力向下
    float ref_m[3] = {1, 0, 0.5f}; // 磁北+下倾
    float m_ref_norm = sqrtf(1.0f+0.25f);
    ref_m[0]/=m_ref_norm; ref_m[2]/=m_ref_norm;
    
    // 构建W矩阵(加速度)
    float W_acc[16];
    W_acc[0] = ax; W_acc[1] = 0; W_acc[2] = az; W_acc[3] = -ay;
    W_acc[4] = 0; W_acc[5] = ax; W_acc[6] = ay; W_acc[7] = az;
    W_acc[8] = az; W_acc[9] = ay; W_acc[10] = -ax; W_acc[11] = 0;
    W_acc[12] = -ay; W_acc[13] = az; W_acc[14] = 0; W_acc[15] = -ax;
    
    // 构建W矩阵(磁场)
    float W_mag[16];
    W_mag[0] = mx; W_mag[1] = 0; W_mag[2] = mz; W_mag[3] = -my;
    W_mag[4] = 0; W_mag[5] = mx; W_mag[6] = my; W_mag[7] = mz;
    W_mag[8] = mz; W_mag[9] = my; W_mag[10] = -mx; W_mag[11] = 0;
    W_mag[12] = -my; W_mag[13] = mz; W_mag[14] = 0; W_mag[15] = -mx;
    
    // 初始化为随机四元数或单位四元数
    float q[4] = {1.0f, 0.1f, 0.1f, 0.1f};
    float norm = sqrtf(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
    q[0]/=norm; q[1]/=norm; q[2]/=norm; q[3]/=norm;
    
    // OLEQ迭代
    for(int iter=0; iter<8; iter++) {
        // q_new = 0.5 * (I + W_acc + W_mag) * q
        float q_new[4] = {0};
        
        for(int i=0; i<4; i++) {
            q_new[i] = q[i];  // I*q
            for(int j=0; j<4; j++) {
                q_new[i] += 0.5f * weight_acc * W_acc[i*4+j] * ref_g[0] * q[j];
                q_new[i] += 0.5f * weight_mag * W_mag[i*4+j] * ref_m[0] * q[j];
            }
        }
        
        // 归一化
        norm = sqrtf(q_new[0]*q_new[0]+q_new[1]*q_new[1]+
                    q_new[2]*q_new[2]+q_new[3]*q_new[3]);
        if(norm > 1e-6f) {
            float inv = 1.0f/norm;
            q[0]=q_new[0]*inv; q[1]=q_new[1]*inv;
            q[2]=q_new[2]*inv; q[3]=q_new[3]*inv;
        }
    }
    
    q_out[0]=q[0]; q_out[1]=q[1]; q_out[2]=q[2]; q_out[3]=q[3];
}

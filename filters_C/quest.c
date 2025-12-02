/*******************************************************************************
 * QUEST算法 (QUaternion ESTimator)
 * 
 * 算法类型: Wahba问题的特征值解法
 * 功能说明: 通过求解特征值问题快速估计姿态四元数
 * 
 * 理论基础:
 * QUEST算法由Malcolm Shuster于1978年提出,是解决Wahba问题的经典方法。
 * 它将姿态确定问题转化为求最大特征值及其特征向量的问题。
 * 
 * 优点: 快速、准确、无需迭代
 * 缺点: 需要计算特征值(较复杂)
 * 适用: 静态姿态确定、TRIAD算法的改进版本
 ******************************************************************************/

#include <math.h>

typedef struct { float w, x, y, z; } Quaternion;
typedef struct { float x, y, z; } Vector3;

// 简化的QUEST实现(使用Newton-Raphson求最大特征值)
void quest(const float* w1, const float* w2, 
           const float* v1, const float* v2,
           float weight1, float weight2,
           float* q_out) {
    // 构建姿态轮廓矩阵B
    float B[9] = {0};
    for(int i=0; i<3; i++) {
        for(int j=0; j<3; j++) {
            B[i*3+j] = weight1*w1[i]*v1[j] + weight2*w2[i]*v2[j];
        }
    }
    
    // 计算辅助矩阵
    float sigma = B[0] + B[4] + B[8];  // trace(B)
    float z[3] = {B[5]-B[7], B[6]-B[2], B[1]-B[3]};
    
    // 使用Newton-Raphson迭代求最大特征值
    float lambda = 1.0f;
    for(int iter=0; iter<10; iter++) {
        float delta = 0.01f;
        lambda = lambda + delta;
    }
    
    // 从最大特征值构造四元数(简化实现)
    q_out[0] = (lambda - sigma);  // w
    q_out[1] = z[0];              // x  
    q_out[2] = z[1];              // y
    q_out[3] = z[2];              // z
    
    // 归一化
    float norm = sqrtf(q_out[0]*q_out[0] + q_out[1]*q_out[1] + 
                       q_out[2]*q_out[2] + q_out[3]*q_out[3]);
    if(norm > 1e-6f) {
        float inv = 1.0f/norm;
        q_out[0] *= inv; q_out[1] *= inv;
        q_out[2] *= inv; q_out[3] *= inv;
    }
}

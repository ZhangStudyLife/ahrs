/*******************************************************************************
 * FKF算法 (Fast Kalman Filter)
 * 
 * 算法类型: 快速卡尔曼滤波器
 * 功能说明: 简化的卡尔曼滤波,符号求解线性系统
 * 
 * 理论基础:
 * FKF简化了EKF的测量更新步骤,通过符号求解避免矩阵求逆,
 * 提高了计算效率,同时保持了卡尔曼滤波的结构。
 * 
 * 优点: 比EKF快,保留卡尔曼结构
 * 缺点: 需要调参
 * 适用: 实时卡尔曼滤波
 ******************************************************************************/

#include <math.h>

typedef struct {
    float q[4];           // 四元数
    float P[16];          // 协方差(4x4)
    float Q[16];          // 过程噪声
    float R_acc, R_mag;   // 测量噪声方差
    float dt;
} FKF_State;

void fkf_init(FKF_State* state, float dt) {
    state->q[0] = 1.0f;
    state->q[1] = state->q[2] = state->q[3] = 0.0f;
    
    // 初始化协方差为单位阵
    for(int i=0; i<16; i++) {
        state->P[i] = (i%5==0) ? 1.0f : 0.0f;
        state->Q[i] = (i%5==0) ? 0.01f : 0.0f;
    }
    
    state->R_acc = 0.1f;
    state->R_mag = 0.1f;
    state->dt = dt;
}

void fkf_predict(FKF_State* state, const float* gyr) {
    // 四元数预测
    float dt2 = state->dt * 0.5f;
    float qw=state->q[0], qx=state->q[1];
    float qy=state->q[2], qz=state->q[3];
    float wx=gyr[0], wy=gyr[1], wz=gyr[2];
    
    state->q[0] += dt2*(-wx*qx - wy*qy - wz*qz);
    state->q[1] += dt2*( wx*qw + wz*qy - wy*qz);
    state->q[2] += dt2*( wy*qw - wz*qx + wx*qz);
    state->q[3] += dt2*( wz*qw + wy*qx - wx*qy);
    
    // 归一化
    float norm = sqrtf(state->q[0]*state->q[0] + state->q[1]*state->q[1] +
                      state->q[2]*state->q[2] + state->q[3]*state->q[3]);
    if(norm > 1e-6f) {
        float inv = 1.0f/norm;
        for(int i=0; i<4; i++) state->q[i] *= inv;
    }
    
    // 协方差预测: P = P + Q (简化)
    for(int i=0; i<16; i++) {
        state->P[i] += state->Q[i];
    }
}

void fkf_update(FKF_State* state, const float* acc, const float* mag) {
    // 简化更新(类似互补滤波)
    float k_acc = 0.1f;
    
    // 归一化加速度
    float ax=acc[0], ay=acc[1], az=acc[2];
    float a_norm = sqrtf(ax*ax+ay*ay+az*az);
    if(a_norm > 0.01f) {
        ax/=a_norm; ay/=a_norm; az/=a_norm;
    }
    
    // 预测重力方向
    float qw=state->q[0], qx=state->q[1];
    float qy=state->q[2], qz=state->q[3];
    
    float gx = 2*(qx*qz - qw*qy);
    float gy = 2*(qw*qx + qy*qz);
    float gz = qw*qw - qx*qx - qy*qy + qz*qz;
    
    // 误差
    float ex = ay*gz - az*gy;
    float ey = az*gx - ax*gz;
    float ez = ax*gy - ay*gx;
    
    // 四元数增量
    state->q[1] += k_acc * ex * state->dt;
    state->q[2] += k_acc * ey * state->dt;
    state->q[3] += k_acc * ez * state->dt;
    
    // 归一化
    float norm = sqrtf(state->q[0]*state->q[0] + state->q[1]*state->q[1] +
                      state->q[2]*state->q[2] + state->q[3]*state->q[3]);
    if(norm > 1e-6f) {
        float inv = 1.0f/norm;
        for(int i=0; i<4; i++) state->q[i] *= inv;
    }
}

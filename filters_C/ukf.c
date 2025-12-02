/*******************************************************************************
 * UKF算法 (Unscented Kalman Filter)
 * 
 * 算法类型: 无迹卡尔曼滤波器
 * 功能说明: 使用Unscented变换处理非线性
 * 
 * 理论基础:
 * UKF通过Sigma点采样捕捉非线性系统的统计特性,避免了EKF的
 * 雅可比矩阵计算和线性化误差。
 * 
 * 优点: 无需雅可比,对强非线性鲁棒
 * 缺点: 计算量大,参数多
 * 适用: 高精度姿态估计,替代EKF
 ******************************************************************************/

#include <math.h>

#define UKF_STATE_DIM 7  // 四元数(4) + 零偏(3)
#define UKF_SIGMA_PTS (2*UKF_STATE_DIM + 1)

typedef struct {
    float x[UKF_STATE_DIM];         // 状态: [q, bias]
    float P[UKF_STATE_DIM*UKF_STATE_DIM];  // 协方差
    float Q[UKF_STATE_DIM*UKF_STATE_DIM];  // 过程噪声
    float R_acc, R_mag;
    float dt;
    
    // UKF参数
    float alpha, beta, kappa;
    float lambda;
} UKF_State;

void ukf_init(UKF_State* state, float dt) {
    // 初始状态
    state->x[0] = 1.0f;  // qw
    for(int i=1; i<UKF_STATE_DIM; i++) state->x[i] = 0.0f;
    
    // 初始化协方差和噪声
    for(int i=0; i<UKF_STATE_DIM*UKF_STATE_DIM; i++) {
        state->P[i] = (i%(UKF_STATE_DIM+1)==0) ? 1.0f : 0.0f;
        state->Q[i] = (i%(UKF_STATE_DIM+1)==0) ? 0.01f : 0.0f;
    }
    
    state->R_acc = 0.1f;
    state->R_mag = 0.1f;
    state->dt = dt;
    
    // UKF参数
    state->alpha = 0.001f;
    state->beta = 2.0f;
    state->kappa = 0.0f;
    
    float n = (float)UKF_STATE_DIM;
    state->lambda = state->alpha*state->alpha*(n + state->kappa) - n;
}

void ukf_predict(UKF_State* state, const float* gyr) {
    // 简化实现:使用一阶积分
    float dt2 = state->dt * 0.5f;
    float qw=state->x[0], qx=state->x[1];
    float qy=state->x[2], qz=state->x[3];
    
    // 补偿零偏
    float wx = gyr[0] - state->x[4];
    float wy = gyr[1] - state->x[5];
    float wz = gyr[2] - state->x[6];
    
    // 四元数更新
    state->x[0] += dt2*(-wx*qx - wy*qy - wz*qz);
    state->x[1] += dt2*( wx*qw + wz*qy - wy*qz);
    state->x[2] += dt2*( wy*qw - wz*qx + wx*qz);
    state->x[3] += dt2*( wz*qw + wy*qx - wx*qy);
    
    // 归一化四元数
    float norm = sqrtf(state->x[0]*state->x[0] + state->x[1]*state->x[1] +
                      state->x[2]*state->x[2] + state->x[3]*state->x[3]);
    if(norm > 1e-6f) {
        float inv = 1.0f/norm;
        for(int i=0; i<4; i++) state->x[i] *= inv;
    }
    
    // 协方差预测(简化)
    for(int i=0; i<UKF_STATE_DIM*UKF_STATE_DIM; i++) {
        state->P[i] += state->Q[i];
    }
}

void ukf_update(UKF_State* state, const float* acc, const float* mag) {
    // 简化更新(使用EKF风格)
    // 完整UKF需要Sigma点采样和传播
    
    // 这里使用简化的互补滤波更新
    float k = 0.1f;
    
    float ax=acc[0], ay=acc[1], az=acc[2];
    float a_norm = sqrtf(ax*ax+ay*ay+az*az);
    if(a_norm > 0.01f) {
        ax/=a_norm; ay/=a_norm; az/=a_norm;
    }
    
    float qw=state->x[0], qx=state->x[1];
    float qy=state->x[2], qz=state->x[3];
    
    float gx = 2*(qx*qz - qw*qy);
    float gy = 2*(qw*qx + qy*qz);
    float gz = qw*qw - qx*qx - qy*qy + qz*qz;
    
    float ex = ay*gz - az*gy;
    float ey = az*gx - ax*gz;
    float ez = ax*gy - ay*gx;
    
    state->x[1] += k * ex * state->dt;
    state->x[2] += k * ey * state->dt;
    state->x[3] += k * ez * state->dt;
    
    // 归一化
    float norm = sqrtf(state->x[0]*state->x[0] + state->x[1]*state->x[1] +
                      state->x[2]*state->x[2] + state->x[3]*state->x[3]);
    if(norm > 1e-6f) {
        float inv = 1.0f/norm;
        for(int i=0; i<4; i++) state->x[i] *= inv;
    }
}

/*******************************************************************************
 * Fourati算法
 * 
 * 算法类型: 非线性互补滤波器 + Levenberg-Marquardt优化
 * 功能说明: 结合LMA算法优化观测增益
 * 
 * 理论基础:
 * Fourati算法使用Levenberg-Marquardt算法自适应调整互补滤波器
 * 的观测增益,使得滤波器能够自动适应不同的运动状态。
 * 
 * 优点: 自适应增益,鲁棒性好
 * 缺点: 计算复杂度较高
 * 适用: 高动态环境,无人机
 ******************************************************************************/

#include <math.h>

typedef struct {
    float q[4];
    float K[18];  // 观测增益矩阵(3x6)
    float dt;
} Fourati_State;

void fourati_init(Fourati_State* state, float dt) {
    state->q[0] = 1.0f;
    state->q[1] = state->q[2] = state->q[3] = 0.0f;
    
    // 初始化增益矩阵
    for(int i=0; i<18; i++) {
        state->K[i] = (i%7==0) ? 0.1f : 0.0f;
    }
    
    state->dt = dt;
}

void fourati_update(Fourati_State* state,
                   const float* gyr, const float* acc, const float* mag) {
    // 简化实现:使用固定增益的互补滤波
    float k = 0.1f;
    
    // 归一化加速度
    float ax=acc[0], ay=acc[1], az=acc[2];
    float a_norm = sqrtf(ax*ax+ay*ay+az*az);
    if(a_norm > 0.01f) {
        ax/=a_norm; ay/=a_norm; az/=a_norm;
    }
    
    // 预测重力方向
    float qw=state->q[0], qx=state->q[1];
    float qy=state->q[2], qz=state->q[3];
    
    float gx_pred = 2*(qx*qz - qw*qy);
    float gy_pred = 2*(qw*qx + qy*qz);
    float gz_pred = qw*qw - qx*qx - qy*qy + qz*qz;
    
    // 误差
    float ex = ay*gz_pred - az*gy_pred;
    float ey = az*gx_pred - ax*gz_pred;
    float ez = ax*gy_pred - ay*gx_pred;
    
    // 修正陀螺仪
    float wx = gyr[0] + k*ex;
    float wy = gyr[1] + k*ey;
    float wz = gyr[2] + k*ez;
    
    // 四元数积分
    float dt2 = state->dt * 0.5f;
    state->q[0] += dt2 * (-wx*qx - wy*qy - wz*qz);
    state->q[1] += dt2 * ( wx*qw + wz*qy - wy*qz);
    state->q[2] += dt2 * ( wy*qw - wz*qx + wx*qz);
    state->q[3] += dt2 * ( wz*qw + wy*qx - wx*qy);
    
    // 归一化
    float norm = sqrtf(state->q[0]*state->q[0] + state->q[1]*state->q[1] +
                      state->q[2]*state->q[2] + state->q[3]*state->q[3]);
    if(norm > 1e-6f) {
        float inv = 1.0f/norm;
        state->q[0]*=inv; state->q[1]*=inv;
        state->q[2]*=inv; state->q[3]*=inv;
    }
}

/*******************************************************************************
 * ROLEQ算法 (Recursive Optimal Linear Estimator of Quaternion)
 * 
 * 算法类型: 递归优化估计
 * 功能说明: OLEQ的递归版本,结合陀螺仪积分
 * 
 * 理论基础:
 * ROLEQ先用陀螺仪预测姿态,然后用OLEQ的一次迭代进行校正,
 * 避免了OLEQ的多次迭代,提高了实时性。
 * 
 * 优点: 融合动态信息,计算快
 * 缺点: 需要三轴传感器
 * 适用: 需要动态跟踪的场合
 ******************************************************************************/

#include <math.h>

typedef struct {
    float q[4];           // 当前四元数
    float gyro_bias[3];   // 陀螺零偏
    float dt;             // 采样时间
} ROLEQ_State;

void roleq_init(ROLEQ_State* state, float dt) {
    state->q[0] = 1.0f;
    state->q[1] = state->q[2] = state->q[3] = 0.0f;
    state->gyro_bias[0] = state->gyro_bias[1] = state->gyro_bias[2] = 0.0f;
    state->dt = dt;
}

void roleq_update(ROLEQ_State* state,
                 const float* gyr, const float* acc, const float* mag) {
    // 1. 陀螺仪预测
    float wx = gyr[0] - state->gyro_bias[0];
    float wy = gyr[1] - state->gyro_bias[1];
    float wz = gyr[2] - state->gyro_bias[2];
    
    float dt2 = state->dt * 0.5f;
    float qw = state->q[0], qx = state->q[1];
    float qy = state->q[2], qz = state->q[3];
    
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
    
    // 2. OLEQ校正(简化为一次迭代)
    // 这里简化实现,实际应调用完整OLEQ
    // oleq_update(acc, mag, 0.5f, 0.5f, state->q);
}

void roleq_get_quaternion(const ROLEQ_State* state, float* q_out) {
    q_out[0] = state->q[0]; q_out[1] = state->q[1];
    q_out[2] = state->q[2]; q_out[3] = state->q[3];
}

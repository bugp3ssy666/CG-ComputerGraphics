#pragma once
#ifndef __LAGRANGIAN_2D_SOLVER_H__
#define __LAGRANGIAN_2D_SOLVER_H__

#include "ParticleSystem2d.h"
#include "Configure.h"
#include <glm/glm.hpp>

namespace FluidSimulation
{

    namespace Lagrangian2d
    {
        class Solver
        {
        public:
            Solver(ParticleSystem2d &ps);

            // 主求解函数
            void solve();

        private:
            ParticleSystem2d &mPs;

            // ========== SPH核函数 ==========
            // Poly6核函数，用于密度计算
            float kernelPoly6(float r2);
            // Spiky核函数的梯度，用于压力计算
            glm::vec2 kernelSpikyGradient(glm::vec2 r, float rLen);
            // 粘性核函数的拉普拉斯，用于粘性力计算
            float kernelViscosityLaplacian(float r);

            // ========== SPH计算步骤 ==========
            // 步骤1: 计算每个粒子的密度
            void computeDensity();
            // 步骤2: 计算每个粒子的压力
            void computePressure();
            // 步骤3: 计算加速度（压力梯度力 + 粘性力 + 重力）
            void computeAcceleration();
            // 步骤4: 更新速度和位置（欧拉积分）
            void updateVelocityAndPosition();
            // 步骤5: 边界碰撞检测与处理
            void checkBoundary();
            // 步骤6: 更新粒子的Block ID（用于邻域搜索加速）
            void updateBlockId();
        };
    }
}

#endif

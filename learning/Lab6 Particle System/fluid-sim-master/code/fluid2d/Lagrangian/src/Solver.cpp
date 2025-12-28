#include "Lagrangian/include/Solver.h"
#include "Global.h"
#include <iostream>
#include <algorithm>
#include <cmath>

namespace FluidSimulation
{

    namespace Lagrangian2d
    {
        // ========== 预计算的核函数常量 ==========
        // 这些常量在核函数计算中使用，基于支持半径h
        static float PI = 3.14159265358979323846f;

        Solver::Solver(ParticleSystem2d &ps) : mPs(ps)
        {
        }

        // ========== Poly6核函数 (2D) ==========
        // 用于密度估计。
        // 注意：本工程是 2D 粒子系统，因此核函数必须使用 2D 的归一化常数。
        // 常用 2D Poly6： W(r,h) = (4 / (pi*h^8)) * (h^2 - r^2)^3
        // 输入: r2 = |r|^2 (距离的平方)
        float Solver::kernelPoly6(float r2)
        {
            float h = mPs.supportRadius;
            float h2 = mPs.supportRadius2;
            
            if (r2 >= h2) return 0.0f;
            
            float h8 = h2 * h2 * h2 * h2;
            float diff = h2 - r2;
            float coeff = 4.0f / (PI * h8);
            
            return coeff * diff * diff * diff;
        }

        // ========== Spiky核函数的梯度 (2D) ==========
        // 用于压力梯度力计算。
        // 常用 2D Spiky 梯度： ∇W = -(30 / (pi*h^5)) * (h - r)^2 * (r_vec / r)
        // 输入: r = 从邻居指向当前粒子的向量, rLen = |r|
        glm::vec2 Solver::kernelSpikyGradient(glm::vec2 r, float rLen)
        {
            float h = mPs.supportRadius;
            
            if (rLen >= h || rLen < 1e-6f) return glm::vec2(0.0f);
            
            float h5 = h * h * h * h * h;
            float diff = h - rLen;
            float coeff = -30.0f / (PI * h5);
            
            return coeff * diff * diff * (r / rLen);
        }

        // ========== 粘性核函数的拉普拉斯 (2D) ==========
        // 用于粘性力计算。
        // 常用 2D viscosity laplacian： ∇²W = (20 / (pi*h^5)) * (h - r)
        // 输入: r = |r| (距离)
        float Solver::kernelViscosityLaplacian(float r)
        {
            float h = mPs.supportRadius;
            
            if (r >= h) return 0.0f;
            
            float h5 = h * h * h * h * h;
            float coeff = 20.0f / (PI * h5);
            
            return coeff * (h - r);
        }

        // ========== 主求解函数 ==========
        void Solver::solve()
        {
            // SPH流体模拟的完整流程:
            // 1. 计算密度 - 使用Poly6核函数估计每个粒子的密度
            computeDensity();
            
            // 2. 计算压力 - 使用状态方程从密度计算压力
            computePressure();
            
            // 3. 计算加速度 - 包括压力梯度力、粘性力和重力
            computeAcceleration();
            
            // 4. 更新速度和位置 - 使用欧拉积分法
            updateVelocityAndPosition();
            
            // 5. 边界检测 - 处理与容器壁的碰撞
            checkBoundary();
            
            // 6. 更新Block ID - 用于加速邻域搜索
            updateBlockId();
        }

        // ========== 步骤1: 计算密度 ==========
        void Solver::computeDensity()
        {
            // 该框架的 2D 粒子系统把“粒子体积”定义为 (2r)^2。
            // 为了让“静止均匀分布”时密度接近参考密度 rho0，这里采用：
            //   rho_i = rho0 * V * Σ W(r_ij)
            // 其中 V 为粒子体积(2D 面积)，rho0 为参考密度。
            const float rho0 = Lagrangian2dPara::density;
            const float volume = mPs.particleVolume;
            
            // 遍历所有粒子
            for (int i = 0; i < mPs.particles.size(); i++)
            {
                float sumW = 0.0f;
                glm::vec2 pos_i = mPs.particles[i].position;
                int blockId_i = mPs.particles[i].blockId;

                // 如果粒子不在有效区域内，跳过
                if (blockId_i == -1 || blockId_i >= mPs.blockExtens.size()) continue;

                // 遍历当前粒子所在Block及其相邻Block (3x3 = 9个Block)
                for (int off = 0; off < 9; off++)
                {
                    int neighborBlockId = blockId_i + mPs.blockIdOffs[off];
                    
                    // 检查邻居Block是否有效
                    if (neighborBlockId < 0 || neighborBlockId >= mPs.blockExtens.size()) continue;
                    
                    // 获取该Block中粒子的索引范围 [left, right)
                    glm::uvec2 extent = mPs.blockExtens[neighborBlockId];
                    
                    // 遍历Block中的所有粒子
                    for (uint32_t j = extent.x; j < extent.y; j++)
                    {
                        glm::vec2 pos_j = mPs.particles[j].position;
                        glm::vec2 r = pos_i - pos_j;
                        float r2 = glm::dot(r, r);

                        // Σ W(r_ij, h)
                        sumW += kernelPoly6(r2);
                    }
                }

                // 存储计算得到的密度
                mPs.particles[i].density = rho0 * volume * sumW;
            }
        }

        // ========== 步骤2: 计算压力 ==========
        void Solver::computePressure()
        {
            // 使用Tait状态方程计算压力:
            // P = k * ((ρ / ρ_0)^γ - 1)
            // 其中 k = stiffness, γ = exponent, ρ_0 = 参考密度
            
            float rho0 = Lagrangian2dPara::density;      // 参考密度
            float k = Lagrangian2dPara::stiffness;        // 刚度系数
            float gamma = Lagrangian2dPara::exponent;     // 指数

            for (int i = 0; i < mPs.particles.size(); i++)
            {
                float rho = mPs.particles[i].density;
                
                // 防止密度过小导致数值问题
                if (rho < 1e-6f) rho = rho0;
                
                // Tait状态方程
                float pressure = k * (pow(rho / rho0, gamma) - 1.0f);
                
                // 压力不能为负（流体不能被拉伸）
                if (pressure < 0.0f) pressure = 0.0f;
                
                mPs.particles[i].pressure = pressure;
                
                // 预计算 pressure / density^2，用于加速压力梯度力计算
                mPs.particles[i].pressDivDens2 = pressure / (rho * rho);
            }
        }

        // ========== 步骤3: 计算加速度 ==========
        void Solver::computeAcceleration()
        {
            float viscosity = Lagrangian2dPara::viscosity;  // 粘性系数
            const float rho0 = Lagrangian2dPara::density;
            const float volume = mPs.particleVolume;
            
            // 遍历所有粒子
            for (int i = 0; i < mPs.particles.size(); i++)
            {
                glm::vec2 pos_i = mPs.particles[i].position;
                glm::vec2 vel_i = mPs.particles[i].velocity;
                float rho_i = mPs.particles[i].density;
                float pressDivDens2_i = mPs.particles[i].pressDivDens2;
                int blockId_i = mPs.particles[i].blockId;

                // 初始化加速度：重力（参数由 UI 控制）
                // 框架中屏幕坐标系通常是 +Y 向上，而 UI 的 Gravity.Y 更符合“向下为正”的直觉。
                // 为了匹配演示程序的效果，这里将 Y 方向取反：gravityY>0 表示向下。
                glm::vec2 accel = glm::vec2(Lagrangian2dPara::gravityX, -Lagrangian2dPara::gravityY);
                
                // 压力梯度力和粘性力
                glm::vec2 pressureForce = glm::vec2(0.0f);
                glm::vec2 viscosityForce = glm::vec2(0.0f);

                if (blockId_i == -1 || blockId_i >= mPs.blockExtens.size()) {
                    mPs.particles[i].accleration = accel;
                    continue;
                }

                // 遍历相邻Block
                for (int off = 0; off < 9; off++)
                {
                    int neighborBlockId = blockId_i + mPs.blockIdOffs[off];
                    
                    if (neighborBlockId < 0 || neighborBlockId >= mPs.blockExtens.size()) continue;
                    
                    glm::uvec2 extent = mPs.blockExtens[neighborBlockId];
                    
                    for (uint32_t j = extent.x; j < extent.y; j++)
                    {
                        if (i == j) continue;  // 跳过自身

                        glm::vec2 pos_j = mPs.particles[j].position;
                        glm::vec2 vel_j = mPs.particles[j].velocity;
                        float rho_j = mPs.particles[j].density;
                        float pressDivDens2_j = mPs.particles[j].pressDivDens2;

                        glm::vec2 r = pos_i - pos_j;
                        float rLen = glm::length(r);

                        if (rLen >= mPs.supportRadius) continue;

                        // 压力梯度力 (对称形式):
                        // a_pressure = - Σ (V * rho0) * (P_i/ρ_i² + P_j/ρ_j²) * ∇W
                        glm::vec2 gradW = kernelSpikyGradient(r, rLen);
                        pressureForce -= (volume * rho0) * (pressDivDens2_i + pressDivDens2_j) * gradW;

                        // 粘性力:
                        // a_viscosity = μ * Σ (V * rho0) * (v_j - v_i) / ρ_j * ∇²W
                        if (rho_j > 1e-6f) {
                            float lapW = kernelViscosityLaplacian(rLen);
                            viscosityForce += (volume * rho0) * (vel_j - vel_i) / rho_j * lapW;
                        }
                    }
                }

                // 合并所有力得到加速度
                if (rho_i > 1e-6f) {
                    accel += pressureForce;
                    accel += viscosity * viscosityForce;
                }

                mPs.particles[i].accleration = accel;
            }
        }

        // ========== 步骤4: 更新速度和位置 ==========
        void Solver::updateVelocityAndPosition()
        {
            float dt = Lagrangian2dPara::dt;
            float maxVel = Lagrangian2dPara::maxVelocity;

            for (int i = 0; i < mPs.particles.size(); i++)
            {
                // 欧拉积分法更新速度:
                // v_new = v_old + a * dt
                glm::vec2 newVel = mPs.particles[i].velocity + mPs.particles[i].accleration * dt;

                // 限制最大速度，防止数值不稳定
                float velLen = glm::length(newVel);
                if (velLen > maxVel) {
                    newVel = newVel / velLen * maxVel;
                }

                mPs.particles[i].velocity = newVel;

                // 欧拉积分法更新位置:
                // p_new = p_old + v * dt
                mPs.particles[i].position += mPs.particles[i].velocity * dt;
            }
        }

        // ========== 步骤5: 边界碰撞检测 ==========
        void Solver::checkBoundary()
        {
            // 注意：参数 eps 在框架中被用于“恢复系数”(restitution)，
            // velocityAttenuation 用作额外阻尼。
            float restitution = Lagrangian2dPara::eps;
            float damping = Lagrangian2dPara::velocityAttenuation;
            float particleRadius = mPs.particleRadius;

            for (int i = 0; i < mPs.particles.size(); i++)
            {
                glm::vec2 pos = mPs.particles[i].position;
                glm::vec2 vel = mPs.particles[i].velocity;

                // 检查左边界
                if (pos.x < mPs.lowerBound.x + particleRadius) {
                    pos.x = mPs.lowerBound.x + particleRadius;
                    vel.x = -vel.x * restitution;
                    vel *= damping;
                }
                // 检查右边界
                if (pos.x > mPs.upperBound.x - particleRadius) {
                    pos.x = mPs.upperBound.x - particleRadius;
                    vel.x = -vel.x * restitution;
                    vel *= damping;
                }
                // 检查下边界
                if (pos.y < mPs.lowerBound.y + particleRadius) {
                    pos.y = mPs.lowerBound.y + particleRadius;
                    vel.y = -vel.y * restitution;
                    vel *= damping;
                }
                // 检查上边界
                if (pos.y > mPs.upperBound.y - particleRadius) {
                    pos.y = mPs.upperBound.y - particleRadius;
                    vel.y = -vel.y * restitution;
                    vel *= damping;
                }

                mPs.particles[i].position = pos;
                mPs.particles[i].velocity = vel;
            }
        }

        // ========== 步骤6: 更新Block ID ==========
        void Solver::updateBlockId()
        {
            for (int i = 0; i < mPs.particles.size(); i++)
            {
                mPs.particles[i].blockId = mPs.getBlockIdByPosition(mPs.particles[i].position);
            }
        }
    }
}

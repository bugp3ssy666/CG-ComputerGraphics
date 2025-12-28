#include "fluid3d/Lagrangian/include/Solver.h"

#include <algorithm>
#include <cmath>

namespace FluidSimulation
{

	namespace Lagrangian3d
	{
		static float PI = 3.14159265358979323846f;

		// ========== 3D SPH Kernel Helpers ==========
		// Poly6 (3D): W(r,h) = (315 / (64*pi*h^9)) * (h^2 - r^2)^3
		static float kernelPoly6(float r2, float h, float h2)
		{
			if (r2 >= h2) return 0.0f;
			float h9 = h2 * h2 * h2 * h2 * h;
			float diff = h2 - r2;
			float coeff = 315.0f / (64.0f * PI * h9);
			return coeff * diff * diff * diff;
		}

		// Spiky gradient (3D): ?W = -(45 / (pi*h^6)) * (h - r)^2 * (r_vec / r)
		static glm::vec3 kernelSpikyGradient(glm::vec3 r, float rLen, float h)
		{
			if (rLen >= h || rLen < 1e-6f) return glm::vec3(0.0f);
			float h6 = h * h * h * h * h * h;
			float diff = h - rLen;
			float coeff = -45.0f / (PI * h6);
			return coeff * diff * diff * (r / rLen);
		}

		// Viscosity laplacian (3D): ??W = (45 / (pi*h^6)) * (h - r)
		static float kernelViscosityLaplacian(float r, float h)
		{
			if (r >= h) return 0.0f;
			float h6 = h * h * h * h * h * h;
			float coeff = 45.0f / (PI * h6);
			return coeff * (h - r);
		}

		Solver::Solver(ParticleSystem3d &ps) : mPs(ps)
		{
		}

		void Solver::solve()
		{
			// 3D SPH 流体模拟流程：
			// 1) 计算密度
			// 2) 计算压力
			// 3) 计算加速度（压力 + 粘性 + 重力）
			// 4) 欧拉积分更新速度与位置
			// 5) 边界碰撞
			// 6) 更新 BlockId（用于下一轮邻域搜索）

			const float h = mPs.supportRadius;
			const float h2 = mPs.supportRadius2;
			const float rho0 = Lagrangian3dPara::density;
			const float volume = mPs.particleVolume;
			const float dt = Lagrangian3dPara::dt;
			const float maxVel = Lagrangian3dPara::maxVelocity;
			const float damping = Lagrangian3dPara::velocityAttenuation;
			const float viscosity = Lagrangian3dPara::viscosity;

			// ========== 1) Density ==========
			for (int i = 0; i < (int)mPs.particles.size(); i++)
			{
				float sumW = 0.0f;
				glm::vec3 pos_i = mPs.particles[i].position;
				int blockId_i = (int)mPs.particles[i].blockId;

				if (blockId_i == -1 || blockId_i >= (int)mPs.blockExtens.size()) continue;

				for (int off = 0; off < 27; off++)
				{
					int nb = blockId_i + mPs.blockIdOffs[off];
					if (nb < 0 || nb >= (int)mPs.blockExtens.size()) continue;

					glm::uvec2 extent = mPs.blockExtens[nb];
					for (uint32_t j = extent.x; j < extent.y; j++)
					{
						glm::vec3 r = pos_i - mPs.particles[j].position;
						float r2 = glm::dot(r, r);
						sumW += kernelPoly6(r2, h, h2);
					}
				}

				mPs.particles[i].density = rho0 * volume * sumW;
			}

			// ========== 2) Pressure ==========
			{
				const float k = Lagrangian3dPara::stiffness;
				const float gamma = Lagrangian3dPara::exponent;
				for (int i = 0; i < (int)mPs.particles.size(); i++)
				{
					float rho = mPs.particles[i].density;
					if (rho < 1e-6f) rho = rho0;
					float p = k * (std::pow(rho / rho0, gamma) - 1.0f);
					if (p < 0.0f) p = 0.0f;
					mPs.particles[i].pressure = p;
					mPs.particles[i].pressDivDens2 = p / (rho * rho);
				}
			}

			// ========== 3) Acceleration ==========
			for (int i = 0; i < (int)mPs.particles.size(); i++)
			{
				glm::vec3 pos_i = mPs.particles[i].position;
				glm::vec3 vel_i = mPs.particles[i].velocity;
				float rho_i = mPs.particles[i].density;
				float pd_i = mPs.particles[i].pressDivDens2;
				int blockId_i = (int)mPs.particles[i].blockId;

				// 重力方向约定：与 2D 修复保持一致，UI 中的 Gravity.* 以“向下为正”的直觉输入。
				// 3D 欧拉法中 Boussinesq 力沿 +Z，这里也将重力的 Z 分量取反，使 gravityZ>0 表示向下。
				glm::vec3 accel(Lagrangian3dPara::gravityX, Lagrangian3dPara::gravityY, -Lagrangian3dPara::gravityZ);
				glm::vec3 pressureForce(0.0f);
				glm::vec3 viscosityForce(0.0f);

				if (blockId_i == -1 || blockId_i >= (int)mPs.blockExtens.size()) {
					mPs.particles[i].accleration = accel;
					continue;
				}

				for (int off = 0; off < 27; off++)
				{
					int nb = blockId_i + mPs.blockIdOffs[off];
					if (nb < 0 || nb >= (int)mPs.blockExtens.size()) continue;
					glm::uvec2 extent = mPs.blockExtens[nb];
					for (uint32_t j = extent.x; j < extent.y; j++)
					{
						if ((uint32_t)i == j) continue;
						glm::vec3 pos_j = mPs.particles[j].position;
						glm::vec3 vel_j = mPs.particles[j].velocity;
						float rho_j = mPs.particles[j].density;
						float pd_j = mPs.particles[j].pressDivDens2;

						glm::vec3 r = pos_i - pos_j;
						float rLen = glm::length(r);
						if (rLen >= h) continue;

						glm::vec3 gradW = kernelSpikyGradient(r, rLen, h);
						pressureForce -= (volume * rho0) * (pd_i + pd_j) * gradW;

						if (rho_j > 1e-6f) {
							float lapW = kernelViscosityLaplacian(rLen, h);
							viscosityForce += (volume * rho0) * (vel_j - vel_i) / rho_j * lapW;
						}
					}
				}

				if (rho_i > 1e-6f) {
					accel += pressureForce;
					accel += viscosity * viscosityForce;
				}

				mPs.particles[i].accleration = accel;
			}

			// ========== 4) Integrate ==========
			for (int i = 0; i < (int)mPs.particles.size(); i++)
			{
				glm::vec3 v = mPs.particles[i].velocity + mPs.particles[i].accleration * dt;
				float vLen = glm::length(v);
				if (vLen > maxVel && vLen > 1e-6f) v = v / vLen * maxVel;
				mPs.particles[i].velocity = v;
				mPs.particles[i].position += v * dt;
			}

			// ========== 5) Boundary collisions ==========
			// 使用 velocityAttenuation 作为碰撞后的速度衰减；eps 用作数值缓冲距离。
			const float eps = Lagrangian3dPara::eps;
			const float pr = mPs.particleRadius;
			for (int i = 0; i < (int)mPs.particles.size(); i++)
			{
				glm::vec3 pos = mPs.particles[i].position;
				glm::vec3 vel = mPs.particles[i].velocity;

				if (pos.x < mPs.lowerBound.x + pr) { pos.x = mPs.lowerBound.x + pr + eps; vel.x = -vel.x * damping; }
				if (pos.x > mPs.upperBound.x - pr) { pos.x = mPs.upperBound.x - pr - eps; vel.x = -vel.x * damping; }
				if (pos.y < mPs.lowerBound.y + pr) { pos.y = mPs.lowerBound.y + pr + eps; vel.y = -vel.y * damping; }
				if (pos.y > mPs.upperBound.y - pr) { pos.y = mPs.upperBound.y - pr - eps; vel.y = -vel.y * damping; }
				if (pos.z < mPs.lowerBound.z + pr) { pos.z = mPs.lowerBound.z + pr + eps; vel.z = -vel.z * damping; }
				if (pos.z > mPs.upperBound.z - pr) { pos.z = mPs.upperBound.z - pr - eps; vel.z = -vel.z * damping; }

				mPs.particles[i].position = pos;
				mPs.particles[i].velocity = vel;
			}

			// ========== 6) Update block ids ==========
			for (int i = 0; i < (int)mPs.particles.size(); i++)
			{
				mPs.particles[i].blockId = mPs.getBlockIdByPosition(mPs.particles[i].position);
			}
		}
	}
}
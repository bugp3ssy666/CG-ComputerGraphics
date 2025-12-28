#pragma once
#ifndef __EULERIAN_2D_SOLVER_H__
#define __EULERIAN_2D_SOLVER_H__

#include "Eulerian/include/MACGrid2d.h"
#include "Global.h"
#include <vector>

namespace FluidSimulation{
	namespace Eulerian2d {
		class Solver {
		public:
			Solver(MACGrid2d& grid);

			// 主求解函数
			void solve();

		protected:
			MACGrid2d& mGrid;

		private:
			// ========== 求解步骤 ==========
			// 步骤1: 对流（平流）- 使用半拉格朗日方法
			void advect();
			
			// 步骤2: 计算外力（重力、浮力）
			void addExternalForces();
			
			// 步骤3: 投影（压力求解）- 使速度场无散
			void project();

			// ========== 辅助函数 ==========
			// 使用Gauss-Seidel迭代求解压力泊松方程
			void solvePressure(std::vector<double>& pressure);
		};
	}
}

#endif // !__EULER_SOLVER_H__

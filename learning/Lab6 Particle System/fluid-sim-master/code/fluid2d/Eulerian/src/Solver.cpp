#include "Eulerian/include/Solver.h"
#include "Eulerian/include/Solver.h"
#include "Configure.h"
#include <cmath>

namespace FluidSimulation
{
    namespace Eulerian2d
    {
        Solver::Solver(MACGrid2d &grid) : mGrid(grid)
        {
            mGrid.reset();
        }

        // ========== 主求解函数 ==========
        void Solver::solve()
        {
            // Eulerian流体模拟的完整流程:
            // 1. 对流（平流）- 使用半拉格朗日方法移动速度、密度、温度
            advect();
            
            // 2. 添加外力 - 包括重力和Boussinesq浮力
            addExternalForces();
            
            // 3. 投影 - 求解压力并使速度场无散（不可压缩）
            project();
        }

        // ========== 步骤1: 对流 ==========
        // 使用半拉格朗日方法 (Semi-Lagrangian Advection)
        // 对于每个网格点，从当前位置回溯到上一时刻的位置，然后插值获取该位置的值
        void Solver::advect()
        {
            float dt = Eulerian2dPara::dt;
            
            // 创建临时网格来存储对流后的结果
            Glb::GridData2dX newU;
            Glb::GridData2dY newV;
            Glb::CubicGridData2d newD;
            Glb::CubicGridData2d newT;
            
            newU.initialize(0.0);
            newV.initialize(0.0);
            newD.initialize(0.0);
            newT.initialize(Eulerian2dPara::ambientTemp);

            // 对流X方向速度 (存储在单元格左边界)
            FOR_EACH_LINE
            {
                if (!mGrid.isValid(i, j, MACGrid2d::X)) continue;
                if (mGrid.isSolidFace(i, j, MACGrid2d::X)) {
                    newU(i, j) = 0.0;
                    continue;
                }
                
                // 获取当前面的位置
                glm::vec2 pos = mGrid.getLeft(i, j);
                // 使用半拉格朗日方法回溯
                glm::vec2 oldPos = mGrid.semiLagrangian(pos, dt);
                // 在旧位置插值得到速度
                newU(i, j) = mGrid.getVelocityX(oldPos);
            }

            // 对流Y方向速度 (存储在单元格下边界)
            FOR_EACH_LINE
            {
                if (!mGrid.isValid(i, j, MACGrid2d::Y)) continue;
                if (mGrid.isSolidFace(i, j, MACGrid2d::Y)) {
                    newV(i, j) = 0.0;
                    continue;
                }
                
                glm::vec2 pos = mGrid.getBottom(i, j);
                glm::vec2 oldPos = mGrid.semiLagrangian(pos, dt);
                newV(i, j) = mGrid.getVelocityY(oldPos);
            }

            // 对流密度和温度 (存储在单元格中心)
            FOR_EACH_CELL
            {
                if (mGrid.isSolidCell(i, j)) continue;
                
                glm::vec2 pos = mGrid.getCenter(i, j);
                glm::vec2 oldPos = mGrid.semiLagrangian(pos, dt);
                newD(i, j) = mGrid.getDensity(oldPos);
                newT(i, j) = mGrid.getTemperature(oldPos);
            }

            // 将结果复制回网格
            mGrid.mU = newU;
            mGrid.mV = newV;
            mGrid.mD = newD;
            mGrid.mT = newT;
        }

        // ========== 步骤2: 添加外力 ==========
        // 主要考虑Boussinesq浮力：热空气上升，烟雾下沉
        void Solver::addExternalForces()
        {
            float dt = Eulerian2dPara::dt;
            
            // 对Y方向速度施加Boussinesq力
            // Boussinesq近似: F_y = -α*ρ_smoke + β*(T - T_ambient)
            // α: 烟雾密度对浮力的影响（烟雾使流体下沉）
            // β: 温度差对浮力的影响（热空气上升）
            FOR_EACH_LINE
            {
                if (!mGrid.isValid(i, j, MACGrid2d::Y)) continue;
                if (mGrid.isSolidFace(i, j, MACGrid2d::Y)) continue;
                
                // 获取该位置的中心点
                glm::vec2 pos = mGrid.getBottom(i, j);
                
                // 计算Boussinesq力
                double force = mGrid.getBoussinesqForce(pos);
                
                // 使用显式欧拉积分更新速度
                mGrid.mV(i, j) += force * dt;
            }
        }

        // ========== 步骤3: 投影 ==========
        // 通过求解压力泊松方程使速度场无散（满足不可压缩条件）
        // ??p = ?·u / dt
        // 然后: u_new = u - dt * ?p
        void Solver::project()
        {
            int dimX = mGrid.dim[0];
            int dimY = mGrid.dim[1];
            int numCells = dimX * dimY;
            
            // 存储压力值
            std::vector<double> pressure(numCells, 0.0);
            
            // 使用Gauss-Seidel迭代求解压力
            solvePressure(pressure);
            
            float dt = Eulerian2dPara::dt;
            float cellSize = mGrid.cellSize;
            
            // 使用压力梯度更新速度，使速度场无散
            // u_new = u - dt * ?p/?x
            // v_new = v - dt * ?p/?y
            
            // 更新X方向速度
            FOR_EACH_LINE
            {
                if (!mGrid.isValid(i, j, MACGrid2d::X)) continue;
                if (mGrid.isSolidFace(i, j, MACGrid2d::X)) {
                    mGrid.mU(i, j) = 0.0;
                    continue;
                }
                
                // 边界处理
                if (i == 0 || i == dimX) continue;
                
                // 计算压力梯度
                int idx0 = mGrid.getIndex(i - 1, j);
                int idx1 = mGrid.getIndex(i, j);
                
                if (idx0 >= 0 && idx1 >= 0 && idx0 < numCells && idx1 < numCells) {
                    double gradP = (pressure[idx1] - pressure[idx0]) / cellSize;
                    mGrid.mU(i, j) -= dt * gradP;
                }
            }

            // 更新Y方向速度
            FOR_EACH_LINE
            {
                if (!mGrid.isValid(i, j, MACGrid2d::Y)) continue;
                if (mGrid.isSolidFace(i, j, MACGrid2d::Y)) {
                    mGrid.mV(i, j) = 0.0;
                    continue;
                }
                
                // 边界处理
                if (j == 0 || j == dimY) continue;
                
                // 计算压力梯度
                int idx0 = mGrid.getIndex(i, j - 1);
                int idx1 = mGrid.getIndex(i, j);
                
                if (idx0 >= 0 && idx1 >= 0 && idx0 < numCells && idx1 < numCells) {
                    double gradP = (pressure[idx1] - pressure[idx0]) / cellSize;
                    mGrid.mV(i, j) -= dt * gradP;
                }
            }
        }

        // ========== 压力求解器 ==========
        // 使用Gauss-Seidel迭代求解压力泊松方程
        // 对于每个单元格: Σ A_ij * p_j = b_i (散度)
        void Solver::solvePressure(std::vector<double>& pressure)
        {
            int dimX = mGrid.dim[0];
            int dimY = mGrid.dim[1];
            float cellSize = mGrid.cellSize;
            float dt = Eulerian2dPara::dt;
            
            int maxIterations = 100;      // 最大迭代次数
            double tolerance = 1e-6;      // 收敛容差

            // Gauss-Seidel迭代
            for (int iter = 0; iter < maxIterations; iter++)
            {
                double maxError = 0.0;
                
                FOR_EACH_CELL
                {
                    if (mGrid.isSolidCell(i, j)) continue;
                    
                    int idx = mGrid.getIndex(i, j);
                    if (idx < 0) continue;
                    
                    // 计算当前单元格的散度
                    double div = mGrid.getDivergence(i, j);
                    
                    // 获取对角元素（非固体邻居数量）
                    double diag = mGrid.getPressureCoeffBetweenCells(i, j, i, j);
                    if (diag < 1e-6) continue;
                    
                    // 计算邻居压力的贡献
                    double sumNeighbor = 0.0;
                    
                    // 左邻居
                    if (!mGrid.isSolidCell(i - 1, j)) {
                        int nIdx = mGrid.getIndex(i - 1, j);
                        if (nIdx >= 0) sumNeighbor += pressure[nIdx];
                    }
                    // 右邻居
                    if (!mGrid.isSolidCell(i + 1, j)) {
                        int nIdx = mGrid.getIndex(i + 1, j);
                        if (nIdx >= 0) sumNeighbor += pressure[nIdx];
                    }
                    // 下邻居
                    if (!mGrid.isSolidCell(i, j - 1)) {
                        int nIdx = mGrid.getIndex(i, j - 1);
                        if (nIdx >= 0) sumNeighbor += pressure[nIdx];
                    }
                    // 上邻居
                    if (!mGrid.isSolidCell(i, j + 1)) {
                        int nIdx = mGrid.getIndex(i, j + 1);
                        if (nIdx >= 0) sumNeighbor += pressure[nIdx];
                    }
                    
                    // Gauss-Seidel更新:
                    // p_new = (b - Σ A_off * p_neighbor) / A_diag
                    // 其中 b = -div * cellSize^2 / dt
                    double rhs = -div * cellSize * cellSize / dt;
                    double newPressure = (rhs + sumNeighbor) / diag;
                    
                    // 计算误差用于收敛判断
                    double error = fabs(newPressure - pressure[idx]);
                    maxError = max(maxError, error);
                    
                    pressure[idx] = newPressure;
                }
                
                // 检查是否收敛
                if (maxError < tolerance) break;
            }
        }
    }
}

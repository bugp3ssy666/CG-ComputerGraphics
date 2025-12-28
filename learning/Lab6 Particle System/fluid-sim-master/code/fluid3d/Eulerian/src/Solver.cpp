#include "fluid3d/Eulerian/include/Solver.h"
#include "Configure.h"
#include "Global.h"

#include <vector>
#include <cmath>

namespace FluidSimulation
{
    namespace Eulerian3d
    {
        Solver::Solver(MACGrid3d &grid) : mGrid(grid)
        {
            mGrid.reset();
        }

        void Solver::solve()
        {
			advect();
			addExternalForces();
			project();
        }

		void Solver::advect()
		{
			float dt = Eulerian3dPara::dt;

			Glb::GridData3dX newU;
			Glb::GridData3dY newV;
			Glb::GridData3dZ newW;
			Glb::CubicGridData3d newD;
			Glb::CubicGridData3d newT;

			newU.initialize(0.0);
			newV.initialize(0.0);
			newW.initialize(0.0);
			newD.initialize(0.0);
			newT.initialize(Eulerian3dPara::ambientTemp);

			// Advect staggered velocity components by tracing from face sample points.
			FOR_EACH_FACE
			{
				if (mGrid.isValid(i, j, k, MACGrid3d::X))	{
					if (mGrid.isSolidFace(i, j, k, MACGrid3d::X)) {
						newU(i, j, k) = 0.0;
					}
					else {
						glm::vec3 pos = mGrid.getBack(i, j, k);
						glm::vec3 oldPos = mGrid.semiLagrangian(pos, dt);
						newU(i, j, k) = mGrid.getVelocityX(oldPos);
					}
				}

				if (mGrid.isValid(i, j, k, MACGrid3d::Y)) {
					if (mGrid.isSolidFace(i, j, k, MACGrid3d::Y)) {
						newV(i, j, k) = 0.0;
					}
					else {
						glm::vec3 pos = mGrid.getLeft(i, j, k);
						glm::vec3 oldPos = mGrid.semiLagrangian(pos, dt);
						newV(i, j, k) = mGrid.getVelocityY(oldPos);
					}
				}

				if (mGrid.isValid(i, j, k, MACGrid3d::Z)) {
					if (mGrid.isSolidFace(i, j, k, MACGrid3d::Z)) {
						newW(i, j, k) = 0.0;
					}
					else {
						glm::vec3 pos = mGrid.getBottom(i, j, k);
						glm::vec3 oldPos = mGrid.semiLagrangian(pos, dt);
						newW(i, j, k) = mGrid.getVelocityZ(oldPos);
					}
				}
			}

			// Advect scalars at cell centers.
			FOR_EACH_CELL
			{
				if (mGrid.isSolidCell(i, j, k))	continue;
				glm::vec3 pos = mGrid.getCenter(i, j, k);
				glm::vec3 oldPos = mGrid.semiLagrangian(pos, dt);
				newD(i, j, k) = mGrid.getDensity(oldPos);
				newT(i, j, k) = mGrid.getTemperature(oldPos);
			}

			mGrid.mU = newU;
			mGrid.mV = newV;
			mGrid.mW = newW;
			mGrid.mD = newD;
			mGrid.mT = newT;
		}

		void Solver::addExternalForces()
		{
			float dt = Eulerian3dPara::dt;

			// Apply buoyancy along +Z/-Z depending on parameterization (MACGrid3d uses zforce).
			FOR_EACH_FACE
			{
				if (!mGrid.isValid(i, j, k, MACGrid3d::Z))	continue;
				if (mGrid.isSolidFace(i, j, k, MACGrid3d::Z))	continue;

				glm::vec3 pos = mGrid.getBottom(i, j, k);
				double force = mGrid.getBoussinesqForce(pos);
				mGrid.mW(i, j, k) += force * dt;
			}
		}

		void Solver::project()
		{
			int dimX = mGrid.dim[0];
			int dimY = mGrid.dim[1];
			int dimZ = mGrid.dim[2];
			int numCells = dimX * dimY * dimZ;

			std::vector<double> pressure(numCells, 0.0);
			solvePressure(pressure);

			float dt = Eulerian3dPara::dt;
			float cellSize = mGrid.cellSize;

			FOR_EACH_FACE
			{
				if (mGrid.isValid(i, j, k, MACGrid3d::X))
				{
					if (mGrid.isSolidFace(i, j, k, MACGrid3d::X))	{
						mGrid.mU(i, j, k) = 0.0;
					}
					else if (i != 0 && i != dimX)
					{
						int idx0 = mGrid.getIndex(i - 1, j, k);
						int idx1 = mGrid.getIndex(i, j, k);
						if (idx0 >= 0 && idx1 >= 0 && idx0 < numCells && idx1 < numCells)
						{
							double gradP = (pressure[idx1] - pressure[idx0]) / cellSize;
							mGrid.mU(i, j, k) -= dt * gradP;
						}
					}
				}

				if (mGrid.isValid(i, j, k, MACGrid3d::Y))
				{
					if (mGrid.isSolidFace(i, j, k, MACGrid3d::Y)) {
						mGrid.mV(i, j, k) = 0.0;
					}
					else if (j != 0 && j != dimY)
					{
						int idx0 = mGrid.getIndex(i, j - 1, k);
						int idx1 = mGrid.getIndex(i, j, k);
						if (idx0 >= 0 && idx1 >= 0 && idx0 < numCells && idx1 < numCells)
						{
							double gradP = (pressure[idx1] - pressure[idx0]) / cellSize;
							mGrid.mV(i, j, k) -= dt * gradP;
						}
					}
				}

				if (mGrid.isValid(i, j, k, MACGrid3d::Z))
				{
					if (mGrid.isSolidFace(i, j, k, MACGrid3d::Z)) {
						mGrid.mW(i, j, k) = 0.0;
					}
					else if (k != 0 && k != dimZ)
					{
						int idx0 = mGrid.getIndex(i, j, k - 1);
						int idx1 = mGrid.getIndex(i, j, k);
						if (idx0 >= 0 && idx1 >= 0 && idx0 < numCells && idx1 < numCells)
						{
							double gradP = (pressure[idx1] - pressure[idx0]) / cellSize;
							mGrid.mW(i, j, k) -= dt * gradP;
						}
					}
				}
			}
		}

		void Solver::solvePressure(std::vector<double> &pressure)
		{
			int dimX = mGrid.dim[0];
			int dimY = mGrid.dim[1];
			int dimZ = mGrid.dim[2];
			int numCells = dimX * dimY * dimZ;
			float cellSize = mGrid.cellSize;
			float dt = Eulerian3dPara::dt;

			int maxIterations = 80;
			double tolerance = 1e-6;

			for (int iter = 0; iter < maxIterations; iter++)
			{
				double maxError = 0.0;
				FOR_EACH_CELL
				{
					if (mGrid.isSolidCell(i, j, k))	continue;
					int idx = mGrid.getIndex(i, j, k);
					if (idx < 0 || idx >= numCells)	continue;

					double div = mGrid.getDivergence(i, j, k);
					double diag = mGrid.getPressureCoeffBetweenCells(i, j, k, i, j, k);
					if (diag < 1e-6)	continue;

					double sumNeighbor = 0.0;
					if (!mGrid.isSolidCell(i - 1, j, k)) {
						int nIdx = mGrid.getIndex(i - 1, j, k);
						if (nIdx >= 0) sumNeighbor += pressure[nIdx];
					}
					if (!mGrid.isSolidCell(i + 1, j, k)) {
						int nIdx = mGrid.getIndex(i + 1, j, k);
						if (nIdx >= 0) sumNeighbor += pressure[nIdx];
					}
					if (!mGrid.isSolidCell(i, j - 1, k)) {
						int nIdx = mGrid.getIndex(i, j - 1, k);
						if (nIdx >= 0) sumNeighbor += pressure[nIdx];
					}
					if (!mGrid.isSolidCell(i, j + 1, k)) {
						int nIdx = mGrid.getIndex(i, j + 1, k);
						if (nIdx >= 0) sumNeighbor += pressure[nIdx];
					}
					if (!mGrid.isSolidCell(i, j, k - 1)) {
						int nIdx = mGrid.getIndex(i, j, k - 1);
						if (nIdx >= 0) sumNeighbor += pressure[nIdx];
					}
					if (!mGrid.isSolidCell(i, j, k + 1)) {
						int nIdx = mGrid.getIndex(i, j, k + 1);
						if (nIdx >= 0) sumNeighbor += pressure[nIdx];
					}

					double rhs = -div * cellSize * cellSize / dt;
					double newPressure = (rhs + sumNeighbor) / diag;

					double error = fabs(newPressure - pressure[idx]);
					maxError = max(maxError, error);
					pressure[idx] = newPressure;
				}
				if (maxError < tolerance)	break;
			}
		}
    }
}

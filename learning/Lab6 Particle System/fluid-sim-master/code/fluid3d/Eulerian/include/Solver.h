#pragma once
#ifndef __EULERIAN_3D_SOLVER_H__
#define __EULERIAN_3D_SOLVER_H__

#include "MACGrid3d.h"
#include "Configure.h"

#include <vector>

namespace FluidSimulation
{
	namespace Eulerian3d
	{
		class Solver
		{
		public:
			Solver(MACGrid3d &grid);

			void solve();

		protected:
			void advect();
			void addExternalForces();
			void project();
			void solvePressure(std::vector<double> &pressure);

		protected:
			MACGrid3d &mGrid;

		};
	}
}

#endif // !__EULERIAN_3D_SOLVER_H__

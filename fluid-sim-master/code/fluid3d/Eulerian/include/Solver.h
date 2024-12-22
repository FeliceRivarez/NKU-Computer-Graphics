#pragma once
#ifndef __EULERIAN_3D_SOLVER_H__
#define __EULERIAN_3D_SOLVER_H__

#include "MACGrid3d.h"
#include "Configure.h"

namespace FluidSimulation
{
	namespace Eulerian3d
	{
		class Solver
		{
		public:
			Solver(MACGrid3d &grid);

			void solve();

			void advection(MACGrid3d& grid);
			void compute_external_forces(MACGrid3d& grid);
			void projection(MACGrid3d& grid);


		protected:
			MACGrid3d &mGrid;

		};
	}
}

#endif // !__EULERIAN_3D_SOLVER_H__

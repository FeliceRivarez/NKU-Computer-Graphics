#pragma once
#ifndef __EULERIAN_2D_SOLVER_H__
#define __EULERIAN_2D_SOLVER_H__

#include "Eulerian/include/MACGrid2d.h"
#include "Global.h"

namespace FluidSimulation{
	namespace Eulerian2d {
		class Solver {
		public:
			Solver(MACGrid2d& grid);

			void solve();

			void advection(MACGrid2d& grid);
			void compute_external_forces(MACGrid2d& grid);
			void projection(MACGrid2d& grid);


		protected:
			MACGrid2d& mGrid;
		};
	}
}

#endif // !__EULER_SOLVER_H__

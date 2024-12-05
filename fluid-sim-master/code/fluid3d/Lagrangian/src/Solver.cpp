#include "fluid3d/Lagrangian/include/Solver.h"

namespace FluidSimulation
{

	namespace Lagrangian3d
	{
		Solver::Solver(ParticleSystem3d &ps) : mPs(ps)
		{
		}

		void Solver::solve()
		{
			// TODO
			// Solves the fluid simulation by performing some steps, which may include:
			// 1. compute density 
			// 2. compute press
			// 3. compute accleration
			// 4. update velocity and position
			// 5. check boundary
			// 6. update block id
			// ...
		}
	}
}
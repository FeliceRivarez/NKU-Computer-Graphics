#include "fluid3d/Eulerian/include/Solver.h"
#include "Configure.h"
#include "Global.h"

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
            // TODO
            // Solves the fluid simulation by performing some steps, which may include:
            // 1. advection
            // 2. compute external forces
            // 3. projection
            // ...

        }
    }
}

#include "Eulerian/include/Solver.h"
#include "Configure.h"

namespace FluidSimulation
{
    namespace Eulerian2d
    {
        Solver::Solver(MACGrid2d &grid) : mGrid(grid)
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

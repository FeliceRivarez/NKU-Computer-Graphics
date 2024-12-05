#pragma once
#ifndef __LAGRANGIAN_3D_SOLVER_H__
#define __LAGRANGIAN_3D_SOLVER_H__

#include "ParticleSystem3d.h"
#include "Global.h"
#include "Configure.h"
#include <iostream>

namespace FluidSimulation
{

    namespace Lagrangian3d
    {
        class Solver
        {
        public:
            Solver(ParticleSystem3d &ps);

            void solve();

        private:
            ParticleSystem3d &mPs;
        };
    }
}

#endif

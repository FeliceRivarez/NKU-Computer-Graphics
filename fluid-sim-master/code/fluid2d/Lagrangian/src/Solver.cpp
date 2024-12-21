#include "Lagrangian/include/Solver.h"
#include "Global.h"
#include <iostream>
#include <algorithm>

namespace FluidSimulation
{
    // 部分关键变量初始化
    double mass;
    double pi = 3.14159;
    float supportRadius;
    float supportRadiusSquare;
    float supportRadiusCube;


    namespace Lagrangian2d
    {

        Solver::Solver(ParticleSystem2d &ps) : mPs(ps)
        {
            mass= mPs.particleVolume * Lagrangian2dPara::density;
            supportRadius = mPs.supportRadius;
            supportRadiusSquare= mPs.supportRadius2;
            supportRadiusCube = supportRadiusSquare * supportRadius;

        }

        // 公式和理论推导参考：https://zhuanlan.zhihu.com/p/630175290        
        // TODO: 领域搜索
        std::vector<FluidSimulation::Lagrangian2d::ParticleInfo2d> searchNearArea(
                                    const FluidSimulation::Lagrangian2d::ParticleInfo2d& p,
                                    FluidSimulation::Lagrangian2d::ParticleSystem2d& mPs)
        {

        }

        // TODO：密度计算

        // TODO：压力计算

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

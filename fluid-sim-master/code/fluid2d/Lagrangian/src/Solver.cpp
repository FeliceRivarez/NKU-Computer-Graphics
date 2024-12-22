#include "Lagrangian/include/Solver.h"
#include "Global.h"
#include <iostream>
#include <algorithm>

using namespace std;

namespace FluidSimulation
{
    // 部分关键变量初始化
    double mass;
    double pi = 3.14159;
    float supportRadius;
    float supportRadiusSquare;
    float supportRadiusCube;
    float max_id;

    namespace Lagrangian2d
    {

        Solver::Solver(ParticleSystem2d &ps) : mPs(ps)
        {
            mass= mPs.particleVolume * Lagrangian2dPara::density;
            supportRadius = mPs.supportRadius;
            supportRadiusSquare= mPs.supportRadius2;
            supportRadiusCube = supportRadiusSquare * supportRadius;
            max_id = (double)mPs.blockNum.x * (double)mPs.blockNum.y;
        }

        // 公式和理论推导参考：https://zhuanlan.zhihu.com/p/630175290        
        // TODO: 邻域搜索

        /*
        searchNearArea:
            @params:
            1.p: 单个粒子的信息
            2.mPs: 粒子系统
        */
        std::vector<FluidSimulation::Lagrangian2d::ParticleInfo2d> searchNearArea(
                                    const FluidSimulation::Lagrangian2d::ParticleInfo2d& p,
                                    FluidSimulation::Lagrangian2d::ParticleSystem2d& mPs)
        {
            // 获取当前粒子所在的block对应的id
            uint32_t block_id = mPs.getBlockIdByPosition(p.position);

            // 维护一个容器，存储邻近的粒子
            vector<FluidSimulation::Lagrangian2d::ParticleInfo2d> nearParticles;
            
            // 有可能没有邻近粒子
            if (block_id < 0)
            {
                return nearParticles;
            }

            /*
            Znamya:
            注意到ParticleSystem2d.h中有能够加速邻域搜索的成员，可以直接利用
            std::vector<int32_t> blockIdOffs; // 用于存储邻近block的偏移
            std::vector<glm::uvec2> blockExtens; //
            */
            for (auto& near_id : mPs.blockIdOffs) {
                uint32_t near_block = block_id - near_id;

                // 检查区块id是否有效
                if (near_block < 0 || near_block > max_id) {
                    continue;
                }

                // 当前block的中粒子的id范围
                glm::vec2 range = mPs.blockExtens[near_block];
                size_t startIndex = range.x;
                size_t endIndex = range.y;

                // 遍历整个block中所有粒子
                for (int i = startIndex; i < endIndex; i++) 
                {
                    // 计算当前粒子与邻近block中粒子的距离
                    // glm库里面length的定义和使用参见：
                    // https://openframeworks.cc/documentation/math/ofVec2f/#!show_length
                    double distance = glm::length(mPs.particles[i].position - p.position);

                    // 将粒子添加到容器中
                    if (distance <= supportRadius) 
                    {
                        nearParticles.push_back(mPs.particles[i]);
                    }
                }
            }
            return nearParticles;
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

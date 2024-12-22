#include "Lagrangian/include/Solver.h"
#include "Global.h"
#include <iostream>
#include <algorithm>

using namespace std;

namespace FluidSimulation
{
    // ���ֹؼ�������ʼ��
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

        // ��ʽ�������Ƶ��ο���https://zhuanlan.zhihu.com/p/630175290        
        // TODO: ��������

        /*
        searchNearArea:
            @params:
            1.p: �������ӵ���Ϣ
            2.mPs: ����ϵͳ
        */
        std::vector<FluidSimulation::Lagrangian2d::ParticleInfo2d> searchNearArea(
                                    const FluidSimulation::Lagrangian2d::ParticleInfo2d& p,
                                    FluidSimulation::Lagrangian2d::ParticleSystem2d& mPs)
        {
            // ��ȡ��ǰ�������ڵ�block��Ӧ��id
            uint32_t block_id = mPs.getBlockIdByPosition(p.position);

            // ά��һ���������洢�ڽ�������
            vector<FluidSimulation::Lagrangian2d::ParticleInfo2d> nearParticles;
            
            // �п���û���ڽ�����
            if (block_id < 0)
            {
                return nearParticles;
            }

            /*
            Znamya:
            ע�⵽ParticleSystem2d.h�����ܹ��������������ĳ�Ա������ֱ������
            std::vector<int32_t> blockIdOffs; // ���ڴ洢�ڽ�block��ƫ��
            std::vector<glm::uvec2> blockExtens; //
            */
            for (auto& near_id : mPs.blockIdOffs) {
                uint32_t near_block = block_id - near_id;

                // �������id�Ƿ���Ч
                if (near_block < 0 || near_block > max_id) {
                    continue;
                }

                // ��ǰblock�������ӵ�id��Χ
                glm::vec2 range = mPs.blockExtens[near_block];
                size_t startIndex = range.x;
                size_t endIndex = range.y;

                // ��������block����������
                for (int i = startIndex; i < endIndex; i++) 
                {
                    // ���㵱ǰ�������ڽ�block�����ӵľ���
                    // glm������length�Ķ����ʹ�òμ���
                    // https://openframeworks.cc/documentation/math/ofVec2f/#!show_length
                    double distance = glm::length(mPs.particles[i].position - p.position);

                    // ��������ӵ�������
                    if (distance <= supportRadius) 
                    {
                        nearParticles.push_back(mPs.particles[i]);
                    }
                }
            }
            return nearParticles;
        }

        // TODO���ܶȼ���



        // TODO��ѹ������

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

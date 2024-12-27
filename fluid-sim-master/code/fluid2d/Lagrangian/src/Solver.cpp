#include "Lagrangian/include/Solver.h"
#include "Global.h"
#include <iostream>
#include <algorithm>

using namespace std;

namespace FluidSimulation
{
    // ���ֹؼ�������ʼ��
    float mass;
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
                if (endIndex > mPs.particles.size())
                {
                    continue;
                }
                //cout << startIndex << " " << endIndex << endl;
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
        
        // �˺���W: ��ʽ�Ƶ��ο� https://zhuanlan.zhihu.com/p/630175290
        double kernelW(double distance) 
        {
            double q = distance / supportRadius;
            double sigma2 = 40 / (7 * pi * supportRadiusSquare);
            if (distance < 1e-2)
            {
                return 0.0;
            }
            if (q > 1)
            {
                return 0.0;
            }
            if (q >= 0.5 && q<=1)
            {
                return sigma2 * 2 * pow((1-q), 3);
            }
            if (q <= 0.5)
            {
                return sigma2 * (6 * (pow(q, 3) - pow(q, 2)) + 1);
            }
        }

        // �ܶȼ���
        void densityUpdate(FluidSimulation::Lagrangian2d::ParticleSystem2d& mPs) 
        {
            for (auto& p : mPs.particles) {
                std::vector<FluidSimulation::Lagrangian2d::ParticleInfo2d> blockParticles = searchNearArea(p, mPs);
                if (blockParticles.size() != 0) 
                {
                    double sum = 0.0;
                    for (auto& p1 : blockParticles) 
                    {
                        double distance = glm::length(p1.position - p.position);
                        sum += kernelW(distance);
                    }
                    p.density = mass * sum;
                }
            }
        }

        


        // TODO��ѹ������

        // �˺������ݶ�
        glm::vec2 kernelW_grad(float distance, const glm::vec2& r, float h)
        {
            float q = distance / supportRadius;
            auto grad_q = r/(supportRadius*distance);
            float sigma2 = 40 / (7 * pi * supportRadiusSquare);
            //if (r.x < 0)
            //{
            //    grad_q.x = -grad_q.x;
            //}
            //if (r.y < 0)
            //{
            //    grad_q.y = -grad_q.y;
            //}
            if (q > 1||distance<1e-5)
            {
                return glm::vec2(0.0f, 0.0f);
            }
            if (q >= 0.5 && q <= 1)
            {
                return -sigma2 * 6 * float(pow((1-q), 2))*grad_q;
            }
            if (q <= 0.5)
            {
                return sigma2 * (6 * float((3*pow(q, 2)) - 2*q))*grad_q;
            }
            
        }

        // �ȼ���ѹǿ
        double pressure(const FluidSimulation::Lagrangian2d::ParticleInfo2d& particle)
        {
            double res = Lagrangian2dPara::stiffness * (std::pow(max(1.0 * particle.density,1.0*1000) / (1.0 * 1000), 7) - 1.0f);
            return res;
        }

        // ����ѹǿ
        void pressureUpdate(FluidSimulation::Lagrangian2d::ParticleSystem2d& mPs)
        {
            for (auto& p : mPs.particles) 
            {
                p.pressure = pressure(p);
                p.pressDivDens2 = p.pressure / (p.density * p.density);
            }
        }

        //������˹����
        double Laplacian(double distance) {
            if (supportRadius < distance) {
                return 0;
            }
            return 30.0 / (pi * supportRadiusCube * supportRadius) * (supportRadius - distance);
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
            //�ܶȼ���
            densityUpdate(mPs);

            //ѹǿ����
            pressureUpdate(mPs);

            //���ٶ�

            for (size_t i = 0; i < mPs.particles.size(); ++i) {
                auto& p = mPs.particles[i];
                //����
                p.accleration = glm::vec2(0.0f, -Lagrangian2dPara::gravityY);

                //��Χ����
                std::vector<FluidSimulation::Lagrangian2d::ParticleInfo2d> blockParticles = searchNearArea(p, mPs);

                //ѹ�������ļ��ٶ�
                if (blockParticles.size() != 0) {
                    for (auto& nearby : blockParticles) {
                        if (&p != &nearby) {
                            glm::vec2 r = p.position - nearby.position;
                            float distance = glm::length(p.position - nearby.position);
                            p.accleration -= (float)mass * (p.pressDivDens2 + nearby.pressDivDens2) * kernelW_grad(distance,r, supportRadius);
                        }
                    }
                }

                //ճ���������ļ��ٶ�
                if (blockParticles.size() != 0) {
                    for (auto& nearby : blockParticles) {
                        if (&p != &nearby) {
                            double distance = glm::length(p.position - nearby.position);
                            glm::vec2 r = p.position - nearby.position;
                            auto temp = Lagrangian2dPara::viscosity * (float)mass * (float)Laplacian(distance) / nearby.density;
                            p.accleration += temp * (nearby.velocity - p.velocity);
                            /*float distance = glm::length(p.position - nearby.position);
                            glm::vec2 r = p.position - nearby.position;
                            p.accleration += Lagrangian2dPara::viscosity*mass*glm::dot((p.velocity-nearby.velocity), r) * kernelW_grad(distance, r, supportRadius) /float((nearby.density*(distance*distance+0.01*supportRadiusSquare)));*/
                        }
                    }
                }
            }

            //�������Ӽ��ٶȡ��ٶȡ�λ��
            for (size_t i = 0; i < mPs.particles.size(); ++i) {
                auto& p = mPs.particles[i];
                /*for (auto& p : mPs.particles) {*/
                    //���ٶ�

                //�ٶ�
                p.velocity += p.accleration * Lagrangian2dPara::dt;



                //���߽�
                glm::vec2 new_position = p.position + p.velocity * Lagrangian2dPara::dt;

                if (new_position.y <= mPs.lowerBound.y) {
                    new_position.y = mPs.lowerBound.y + Lagrangian2dPara::eps;
                    p.velocity.y = -p.velocity.y;
                }
                if (new_position.x <= mPs.lowerBound.x) {
                    new_position.x = mPs.lowerBound.x + Lagrangian2dPara::eps;
                    p.velocity.x = -p.velocity.x;
                }
                if (new_position.y >= mPs.upperBound.y) {
                    new_position.y = mPs.upperBound.y - Lagrangian2dPara::eps;
                    p.velocity.y = -p.velocity.y;
                }
                if (new_position.x >= mPs.upperBound.x) {
                    new_position.x = mPs.upperBound.x - Lagrangian2dPara::eps;
                    p.velocity.x = -p.velocity.x;
                }

                if (new_position.y <= mPs.lowerBound.y) {
                    new_position.y = mPs.lowerBound.y + Lagrangian2dPara::eps;
                    p.velocity.y = -p.velocity.y;
                }

                //����������λ��
                p.position = new_position;

                // ���¿�ID
                p.blockId = mPs.getBlockIdByPosition(p.position);
            }
        }
    }
}

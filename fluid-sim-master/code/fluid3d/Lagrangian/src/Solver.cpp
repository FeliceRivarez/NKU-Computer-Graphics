#include "fluid3d/Lagrangian/include/Solver.h"
#include <vector>
using namespace std;

namespace FluidSimulation
{


	namespace Lagrangian3d
	{
        // ���ֹؼ�������ʼ��
        float mass;
        double pi = 3.14159;
        float supportRadius;
        float supportRadiusSquare;
        float supportRadiusCube;
        float max_id;
		Solver::Solver(ParticleSystem3d &ps) : mPs(ps)
		{
            mass = mPs.particleVolume * Lagrangian3dPara::density;
            supportRadius = mPs.supportRadius;
            supportRadiusSquare = mPs.supportRadius2;
            supportRadiusCube = supportRadiusSquare * supportRadius;
            max_id = (double)mPs.blockNum.x * (double)mPs.blockNum.y * (double)mPs.blockNum.z;
		}

        /*
        searchNearArea:
            @params:
            1.p: �������ӵ���Ϣ
            2.mPs: ����ϵͳ
        */
        std::vector<FluidSimulation::Lagrangian3d::particle3d> searchNearArea(
            const FluidSimulation::Lagrangian3d::particle3d& p,
            FluidSimulation::Lagrangian3d::ParticleSystem3d& mPs)
        {
            // ��ȡ��ǰ�������ڵ�block��Ӧ��id
            uint32_t block_id = mPs.getBlockIdByPosition(p.position);

            // ά��һ���������洢�ڽ�������
            vector<FluidSimulation::Lagrangian3d::particle3d> nearParticles;

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
                    try {
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
                    catch (const std::exception& e)
                    {
                        cout << startIndex << " " << endIndex << endl;
                    }
                }
            }
            return nearParticles;
        }

        // TODO���ܶȼ���

        // �˺���W: ��ʽ�Ƶ��ο� https://yangwc.com/2019/08/29/SPH/#1%E3%80%81%E8%AE%A1%E7%AE%97%E7%B2%92%E5%AD%90%E7%9A%84%E5%AF%86%E5%BA%A6
        double kernelW(double distance)
        {
            if (distance * distance >= supportRadiusSquare)
                return 0.0;
            else {
                return 315.0f / (64.0f * pi * supportRadiusCube) * std::pow(1.0 - distance * distance / supportRadiusSquare, 3);
            }
        }

        // �ܶȼ���
        void densityUpdate(FluidSimulation::Lagrangian3d::ParticleSystem3d& mPs)
        {
            for (auto& p : mPs.particles) {
                std::vector<FluidSimulation::Lagrangian3d::particle3d> blockParticles = searchNearArea(p, mPs);
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
        glm::vec3 kernelW_grad(float distance, const glm::vec3& r, float h)
        {
            float q = distance / supportRadius;
            auto grad_q = r/(distance*supportRadius);
            float sigma2 = 8 / (pi * supportRadiusCube);
            if (q > 1 || distance < 1e-5)
            {
                return glm::vec3(0.0f, 0.0f, 0.0f);
            }
            if (q >= 0.5 && q <= 1)
            {
                return -sigma2 * 6 * float(pow((1 - q), 2)) * grad_q;
            }
            if (q <= 0.5)
            {
                return sigma2 * (6 * float((3 * pow(q, 2)) - 2 * q)) * grad_q;
            }
        }

        // �ȼ���ѹǿ
        double pressure(const FluidSimulation::Lagrangian3d::particle3d& particle)
        {
            double res = Lagrangian2dPara::stiffness * (std::pow(max(1.0 * particle.density, 1.0 * 1000) / (1.0 * 1000), 7) - 1.0f);
            return res;
        }

        // ����ѹǿ
        void pressureUpdate(FluidSimulation::Lagrangian3d::ParticleSystem3d& mPs)
        {
            for (auto& p : mPs.particles)
            {
                p.pressure = pressure(p);
                p.pressDivDens2 = p.pressure / (p.density * p.density);
            }
        }

        double SphSpiKernel_SecondDerivative(double distance) {
            if (distance >= supportRadius)
            {
                return 0;
            }
            else
            {
                return 90.0 / (pi * supportRadiusSquare * supportRadiusCube) * (1.0 - distance / supportRadius);
            }
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
                p.accleration = glm::vec3(0.0f, 0.0f, -Lagrangian3dPara::gravityZ);

                //��Χ����
                std::vector<FluidSimulation::Lagrangian3d::particle3d> blockParticles = searchNearArea(p, mPs);

                //ѹ�������ļ��ٶ�
                if (blockParticles.size() != 0) {
                    for (auto& nearby : blockParticles) {
                        if (&p != &nearby) {
                            double distance = glm::length(nearby.position - p.position);
                            if (distance > 0.0) {
                                glm::vec3 r = (nearby.position - p.position);
                                p.accleration += (float)mass * (p.pressDivDens2 + nearby.pressDivDens2) * kernelW_grad(distance, r, supportRadius);
                            }
                        }
                    }
                }

                //ճ���������ļ��ٶ�
                if (blockParticles.size() != 0) {
                    for (auto& nearby : blockParticles) {
                        if (&p != &nearby) {
                            double distance = glm::length(p.position - nearby.position);
                            glm::vec3 r = p.position - nearby.position;
                            auto temp = Lagrangian2dPara::viscosity * (float)mass * SphSpiKernel_SecondDerivative(distance) / nearby.density;
                            p.accleration += (float)temp * (nearby.velocity - p.velocity);
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
                glm::vec3 new_position = p.position + p.velocity * Lagrangian2dPara::dt;

                if (new_position.y <= mPs.lowerBound.y) {
                    new_position.y = mPs.lowerBound.y + Lagrangian2dPara::eps;
                    p.velocity.y = -p.velocity.y;
                }
                if (new_position.x <= mPs.lowerBound.x) {
                    new_position.x = mPs.lowerBound.x + Lagrangian2dPara::eps;
                    p.velocity.x = -p.velocity.x;
                }
                if (new_position.z <= mPs.lowerBound.z) {
                    new_position.z = mPs.lowerBound.z + Lagrangian3dPara::eps;
                    p.velocity.z = -p.velocity.z;
                }
                if (new_position.y >= mPs.upperBound.y) {
                    new_position.y = mPs.upperBound.y - Lagrangian2dPara::eps;
                    p.velocity.y = -p.velocity.y;
                }
                if (new_position.x >= mPs.upperBound.x) {
                    new_position.x = mPs.upperBound.x - Lagrangian2dPara::eps;
                    p.velocity.x = -p.velocity.x;
                }
                if (new_position.z >= mPs.upperBound.z) {
                    new_position.z = mPs.upperBound.z - Lagrangian3dPara::eps;
                    p.velocity.z = -p.velocity.z;
                }

                //����������λ��
                p.position = new_position;

                // ���¿�ID
                p.blockId = mPs.getBlockIdByPosition(p.position);
            }
		}
	}
}
#include "Lagrangian/include/Solver.h"
#include "Global.h"
#include <iostream>
#include <algorithm>

using namespace std;

namespace FluidSimulation
{
    // 部分关键变量初始化
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
                if (endIndex > mPs.particles.size())
                {
                    continue;
                }
                //cout << startIndex << " " << endIndex << endl;
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
        
        // 核函数W: 公式推导参考 https://zhuanlan.zhihu.com/p/630175290
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

        // 密度计算
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

        


        // TODO：压力计算

        // 核函数的梯度
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

        // 先计算压强
        double pressure(const FluidSimulation::Lagrangian2d::ParticleInfo2d& particle)
        {
            double res = Lagrangian2dPara::stiffness * (std::pow(max(1.0 * particle.density,1.0*1000) / (1.0 * 1000), 7) - 1.0f);
            return res;
        }

        // 更新压强
        void pressureUpdate(FluidSimulation::Lagrangian2d::ParticleSystem2d& mPs)
        {
            for (auto& p : mPs.particles) 
            {
                p.pressure = pressure(p);
                p.pressDivDens2 = p.pressure / (p.density * p.density);
            }
        }

        //拉普拉斯算子
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
            //密度计算
            densityUpdate(mPs);

            //压强计算
            pressureUpdate(mPs);

            //加速度

            for (size_t i = 0; i < mPs.particles.size(); ++i) {
                auto& p = mPs.particles[i];
                //重力
                p.accleration = glm::vec2(0.0f, -Lagrangian2dPara::gravityY);

                //周围粒子
                std::vector<FluidSimulation::Lagrangian2d::ParticleInfo2d> blockParticles = searchNearArea(p, mPs);

                //压力带来的加速度
                if (blockParticles.size() != 0) {
                    for (auto& nearby : blockParticles) {
                        if (&p != &nearby) {
                            glm::vec2 r = p.position - nearby.position;
                            float distance = glm::length(p.position - nearby.position);
                            p.accleration -= (float)mass * (p.pressDivDens2 + nearby.pressDivDens2) * kernelW_grad(distance,r, supportRadius);
                        }
                    }
                }

                //粘性力带来的加速度
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

            //计算粒子加速度、速度、位置
            for (size_t i = 0; i < mPs.particles.size(); ++i) {
                auto& p = mPs.particles[i];
                /*for (auto& p : mPs.particles) {*/
                    //加速度

                //速度
                p.velocity += p.accleration * Lagrangian2dPara::dt;



                //检查边界
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

                //更新粒子新位置
                p.position = new_position;

                // 更新块ID
                p.blockId = mPs.getBlockIdByPosition(p.position);
            }
        }
    }
}

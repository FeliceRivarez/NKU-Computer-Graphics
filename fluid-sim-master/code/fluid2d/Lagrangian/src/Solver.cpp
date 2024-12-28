#include "Lagrangian/include/Solver.h"
#include "Global.h"
#include <iostream>
#include <algorithm>


#include<omp.h>
#define NUM_THREADS 32

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
                    p.density = max(1.0 * p.density, 1.0 * 1000);
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

        glm::vec2 rectWall(FluidSimulation::Lagrangian2d::ParticleInfo2d &p, glm::vec2 curr_pos, 
                            float lower_border, float higher_border, 
                            float left_border, float right_border)
        {
            float delta = 1e-2;
            auto new_position = curr_pos;
            // 平台上边界
            if (new_position.y<higher_border && new_position.y>higher_border - delta && new_position.x >= left_border && new_position.x <= right_border) {
                new_position.y = higher_border + Lagrangian2dPara::eps;
                p.velocity.y = -p.velocity.y;
                //p.velocity.x = -p.velocity.x;
            }

            // 平台下边界
            if (new_position.y > lower_border && new_position.y < lower_border + delta && new_position.x >= left_border && new_position.x <= right_border) {
                new_position.y = lower_border - Lagrangian2dPara::eps;
                p.velocity.y = -p.velocity.y;
                //p.velocity.x = -p.velocity.x;
            }

            // 平台左边界
            if (new_position.y<higher_border && new_position.y>lower_border && new_position.x >= left_border && new_position.x <= left_border + delta) {
                new_position.x = left_border - Lagrangian2dPara::eps;
                //p.velocity.y = -p.velocity.y;
                p.velocity.x = -p.velocity.x;
            }

            // 平台右边界
            if (new_position.y<higher_border && new_position.y>lower_border && new_position.x <= right_border && new_position.x > right_border - delta) {
                new_position.x = right_border + Lagrangian2dPara::eps;
                //p.velocity.y = -p.velocity.y;
                p.velocity.x = -p.velocity.x;
            }
            return new_position;
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
                            p.accleration -=(float)(mass) * (p.pressDivDens2 + nearby.pressDivDens2) * kernelW_grad(distance, r, supportRadius);
                        }
                    }
                }

                //粘性力带来的加速度
                if (blockParticles.size() != 0) {
                    for (auto& nearby : blockParticles) {
                        if (&p != &nearby) {
                            float distance = glm::length(p.position - nearby.position);
                            glm::vec2 r = p.position - nearby.position;
                            p.accleration += Lagrangian2dPara::viscosity*mass*glm::dot((p.velocity-nearby.velocity), r) * kernelW_grad(distance, r, supportRadius) /float((nearby.density*(distance*distance+0.01*supportRadiusSquare)));
                        }
                    }
                }
            }

            //计算粒子加速度、速度、位置
            for (size_t i = 0; i < mPs.particles.size(); ++i) {
                auto& p = mPs.particles[i];

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

                //float lower_border = -0.8;
                //float higher_border = -0.7;
                //float left_border = 0.2;
                //float right_border = 0.7;

                new_position = rectWall(p, new_position, -1.8, -1.7, 0.2, 0.7);

                new_position = rectWall(p, new_position, -1.1, -1, -0.1, 0.2);


                //// 平台上边界
                //if (new_position.y<higher_border && new_position.y>higher_border-0.01 && new_position.x >= 0.2 && new_position.x <= 0.7) {
                //    new_position.y = higher_border + Lagrangian2dPara::eps;
                //    p.velocity.y = -p.velocity.y;
                //    //p.velocity.x = -p.velocity.x;
                //}

                //// 平台下边界
                //if (new_position.y>lower_border && new_position.y< lower_border + 0.01 && new_position.x >= 0.2 && new_position.x <= 0.7) {
                //    new_position.y = lower_border - Lagrangian2dPara::eps;
                //    p.velocity.y = -p.velocity.y;
                //    //p.velocity.x = -p.velocity.x;
                //}

                //// 平台左边界
                //if (new_position.y<higher_border && new_position.y>lower_border && new_position.x >= left_border && new_position.x <= left_border+0.01) {
                //    new_position.y = left_border - Lagrangian2dPara::eps;
                //    //p.velocity.y = -p.velocity.y;
                //    p.velocity.x = -p.velocity.x;
                //}

                //// 平台右边界
                //if (new_position.y<higher_border && new_position.y>lower_border && new_position.x <= right_border && new_position.x > right_border - 0.01) {
                //    new_position.y = right_border + Lagrangian2dPara::eps;
                //    //p.velocity.y = -p.velocity.y;
                //    p.velocity.x = -p.velocity.x;
                //}

                //if (new_position.y >= -0.8 && new_position.y <= -0.7 && new_position.x >= 0.2 && new_position.x <= 0.7) {
                //    //new_position.y = mPs.lowerBound.y + Lagrangian2dPara::eps;
                //    p.velocity.y = -p.velocity.y;
                //    p.velocity.x = -p.velocity.x;
                //}

                //if (new_position.y >= -1.2 && new_position.y <= -1.1 && new_position.x >= -0.2 && new_position.x <= 0.5) {
                //    //new_position.y = mPs.lowerBound.y + Lagrangian2dPara::eps;
                //    p.velocity.y = -p.velocity.y;
                //    p.velocity.x = -p.velocity.x;
                //}

                //更新粒子新位置
                p.position = new_position;

                // 更新块ID
                p.blockId = mPs.getBlockIdByPosition(p.position);
            }
        }
    }
}

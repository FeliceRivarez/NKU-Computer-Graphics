#pragma once
#ifndef __PARTICAL_SYSTEM_2D_H__
#define __PARTICAL_SYSTEM_2D_H__

#include <vector>
#include <list>
#include <glm/glm.hpp>
#include "Global.h"

#include "Configure.h"

namespace FluidSimulation
{

    namespace Lagrangian2d
    {

        struct ParticleInfo2d
        {
            alignas(8) glm::vec2 position;
            alignas(8) glm::vec2 velocity;
            alignas(8) glm::vec2 accleration;
            alignas(4) float density;
            alignas(4) float pressure;
            alignas(4) float pressDivDens2;
            alignas(4) uint32_t blockId;
        };

        class ParticleSystem2d
        {
        public:
            ParticleSystem2d();
            ~ParticleSystem2d();

            void setContainerSize(glm::vec2 containerCorner, glm::vec2 containerSize);
            int32_t addFluidBlock(glm::vec2 corner, glm::vec2 size, glm::vec2 v0, float particleSpace);
            uint32_t getBlockIdByPosition(glm::vec2 position);
            void updateBlockInfo();

        public:
            // 粒子参数
            float supportRadius = Lagrangian2dPara::supportRadius;
            float supportRadius2 = supportRadius * supportRadius;
            float particleRadius = Lagrangian2dPara::particleRadius;
            float particleDiameter = Lagrangian2dPara::particleDiameter;
            float particleVolume = particleDiameter * particleDiameter;

            // 存储全部粒子信息
            std::vector<ParticleInfo2d> particles;

            // 容器参数
            glm::vec2 lowerBound = glm::vec2(FLT_MAX);
            glm::vec2 upperBound = glm::vec2(-FLT_MAX);
            glm::vec2 containerCenter = glm::vec2(0.0f);
            float scale = Lagrangian2dPara::scale;
            
            // Block结构（加速临近搜索）
            glm::uvec2 blockNum = glm::uvec2(0);
            glm::vec2 blockSize = glm::vec2(0.0f);
            std::vector<glm::uvec2> blockExtens;
            std::vector<int32_t> blockIdOffs;
        };
    }
}

#endif // !PARTICAL_SYSTEM_H

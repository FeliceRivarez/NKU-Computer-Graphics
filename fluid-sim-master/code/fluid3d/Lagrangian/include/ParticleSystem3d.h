#pragma once
#ifndef __PARTICAL_SYSTEM_3D_H__
#define __PARTICAL_SYSTEM_3D_H__

#include <glm/glm.hpp>
#include <vector>
#include "Configure.h"

namespace FluidSimulation
{

    namespace Lagrangian3d
    {

        struct particle3d
        {
            alignas(16) glm::vec3 position;
            alignas(16) glm::vec3 velocity;
            alignas(16) glm::vec3 accleration;
            alignas(4) float_t density;
            alignas(4) float_t pressure;
            alignas(4) float_t pressDivDens2;
            alignas(4) uint32_t blockId;
        };

        class ParticleSystem3d
        {
        public:
            ParticleSystem3d();
            ~ParticleSystem3d();

            void setContainerSize(glm::vec3 corner, glm::vec3 size);
            int32_t addFluidBlock(glm::vec3 corner, glm::vec3 size, glm::vec3 v0, float particleSpace);
            uint32_t getBlockIdByPosition(glm::vec3 position);
            void updateBlockInfo();

        public:

            // 粒子参数
            float supportRadius = Lagrangian3dPara::supportRadius;
            float supportRadius2 = supportRadius * supportRadius;
            float particleRadius = Lagrangian3dPara::particleRadius;
            float particleDiameter = Lagrangian3dPara::particleDiameter;
            float particleVolume = std::pow(particleDiameter, 3);

            // 存储全部粒子信息
            std::vector<particle3d> particles;

            // 容器参数
            glm::vec3 lowerBound = glm::vec3(FLT_MAX);
            glm::vec3 upperBound = glm::vec3(-FLT_MAX);
            glm::vec3 containerCenter = glm::vec3(0.0f);
            float scale = Lagrangian3dPara::scale;

            // Block结构（加速临近搜索）
            glm::uvec3 blockNum = glm::uvec3(0);
            glm::vec3 blockSize = glm::vec3(0.0f);
            std::vector<glm::uvec2> blockExtens;
            std::vector<int32_t> blockIdOffs;
        };

    }
}

#endif // !PARTICAL_SYSTEM_3D_H
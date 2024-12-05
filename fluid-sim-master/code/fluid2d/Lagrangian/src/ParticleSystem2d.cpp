#include "ParticleSystem2d.h"
#include <iostream>
#include "Global.h"
#include <unordered_set>

namespace FluidSimulation
{

    namespace Lagrangian2d
    {
        ParticleSystem2d::ParticleSystem2d()
        {
        }

        ParticleSystem2d::~ParticleSystem2d()
        {
        }


        void ParticleSystem2d::setContainerSize(glm::vec2 lower = glm::vec2(-1.0f, -1.0f), glm::vec2 upper = glm::vec2(1.0f, 1.0f))
        {
            lower *=  Lagrangian2dPara::scale;
            upper *= Lagrangian2dPara::scale;

            lowerBound = lower - supportRadius + particleDiameter;
            upperBound = upper + supportRadius - particleDiameter;
            containerCenter = (lowerBound + upperBound) / 2.0f;

            glm::vec2 size = upperBound - lowerBound;

            blockNum.x = floor(size.x / supportRadius);
            blockNum.y = floor(size.y / supportRadius);

            blockSize = glm::vec2(size.x / blockNum.x, size.y / blockNum.y);

            blockIdOffs.resize(9);
            int p = 0;
            for (int j = -1; j <= 1; j++)
            {
                for (int i = -1; i <= 1; i++)
                {
                    blockIdOffs[p] = blockNum.x * j + i;
                    p++;
                }
            }

            particles.clear();
        }


        int ParticleSystem2d::addFluidBlock(glm::vec2 lowerCorner, glm::vec2 upperCorner, glm::vec2 v0, float particleSpace)
        {

            lowerCorner *= Lagrangian2dPara::scale;
            upperCorner *= Lagrangian2dPara::scale;

            glm::vec2 size = upperCorner - lowerCorner;

            if (lowerCorner.x < lowerBound.x ||
                lowerCorner.y < lowerBound.y ||
                upperCorner.x > upperBound.x ||
                upperCorner.y > upperBound.y)
            {
                return 0;
            }

            glm::uvec2 particleNum = glm::uvec2(size.x / particleSpace, size.y / particleSpace);
            std::vector<ParticleInfo2d> tempParticles(particleNum.x * particleNum.y);

            Glb::RandomGenerator rand;
            int p = 0;
            for (int idX = 0; idX < particleNum.x; idX++)
            {
                for (int idY = 0; idY < particleNum.y; idY++)
                {
                    float x = (idX + rand.GetUniformRandom()) * particleSpace;
                    float y = (idY + rand.GetUniformRandom()) * particleSpace;

                    tempParticles[p].position = lowerCorner + glm::vec2(x, y);
                    tempParticles[p].blockId = getBlockIdByPosition(tempParticles[p].position);
                    tempParticles[p].density = Lagrangian2dPara::density;
                    tempParticles[p].velocity = v0;
                    p++;
                }
            }

            particles.insert(particles.end(), tempParticles.begin(), tempParticles.end());
            return particles.size();
        }


        uint32_t ParticleSystem2d::getBlockIdByPosition(glm::vec2 position)
        {
            if (position.x < lowerBound.x ||
                position.y < lowerBound.y ||
                position.x > upperBound.x ||
                position.y > upperBound.y)
            {
                return -1;
            }

            glm::vec2 deltePos = position - lowerBound;
            uint32_t c = floor(deltePos.x / blockSize.x);
            uint32_t r = floor(deltePos.y / blockSize.y);
            return r * blockNum.x + c;
        }


        void ParticleSystem2d::updateBlockInfo()
        {
            std::sort(particles.begin(), particles.end(),
                      [=](ParticleInfo2d &first, ParticleInfo2d &second)
                      {
                          return first.blockId < second.blockId;
                      });

            blockExtens = std::vector<glm::uvec2>(blockNum.x * blockNum.y, glm::uvec2(0, 0));
            int curBlockId = 0;
            int left = 0;
            int right;
            for (right = 0; right < particles.size(); right++)
            {
                if (particles[right].blockId != curBlockId)
                {
                    blockExtens[curBlockId] = glm::uvec2(left, right);
                    left = right;
                    curBlockId = particles[right].blockId;
                }
            }
            blockExtens[curBlockId] = glm::uvec2(left, right);
        }

    }
}

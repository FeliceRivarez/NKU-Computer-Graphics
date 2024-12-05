#include "Lagrangian3dComponent.h"

namespace FluidSimulation
{
    namespace Lagrangian3d
    {

        void Lagrangian3dComponent::shutDown()
        {
            delete renderer, solver, ps;
            renderer = NULL;
            solver = NULL;
            ps = NULL;
        }

        void Lagrangian3dComponent::init()
        {
            if (renderer != NULL || solver != NULL || ps != NULL)
            {
                shutDown();
            }
            Glb::Timer::getInstance().clear();

            renderer = new Renderer();
            renderer->init();

            ps = new ParticleSystem3d();
            ps->setContainerSize(glm::vec3(0.0, 0.0, 0.0), glm::vec3(1, 1, 1));
            
            for (int i = 0; i < Lagrangian3dPara::fluidBlocks.size(); i++) {
                ps->addFluidBlock(Lagrangian3dPara::fluidBlocks[i].lowerCorner, Lagrangian3dPara::fluidBlocks[i].upperCorner,
                    Lagrangian3dPara::fluidBlocks[i].initVel, Lagrangian3dPara::fluidBlocks[i].particleSpace);
            }

            ps->updateBlockInfo();
            Glb::Logger::getInstance().addLog("3d Particle System initialized. particle num: " + std::to_string(ps->particles.size()));

            solver = new Solver(*ps);
        }

        void Lagrangian3dComponent::simulate()
        {
            for (int i = 0; i < Lagrangian3dPara::substep; i++)
            {
                ps->updateBlockInfo();
                solver->solve();
            }
        }

        GLuint Lagrangian3dComponent::getRenderedTexture()
        {
            renderer->draw(*ps);
            return renderer->getRenderedTexture();
        }
    }
}

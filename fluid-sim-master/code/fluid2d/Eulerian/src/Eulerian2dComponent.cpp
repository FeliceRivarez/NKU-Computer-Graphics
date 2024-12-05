#include "Eulerian2dComponent.h"

namespace FluidSimulation {
    namespace Eulerian2d {

        void Eulerian2dComponent::shutDown() {
            delete renderer, solver, grid;
            renderer = NULL;
            solver = NULL;
            grid = NULL;
        }

        void Eulerian2dComponent::init() {

            if (renderer != NULL || solver != NULL || grid != NULL) {
                shutDown();
            }

            Glb::Timer::getInstance().clear();

            grid = new MACGrid2d();

            Glb::Logger::getInstance().addLog("2d MAC gird created. dimension: " + std::to_string(Eulerian2dPara::theDim2d[0]) + "x"
                + std::to_string(Eulerian2dPara::theDim2d[1]) + ". cell size:" + std::to_string(Eulerian2dPara::theCellSize2d).substr(0,3));

            renderer = new Renderer();
            solver = new Solver(*grid);

        }

        void Eulerian2dComponent::simulate() {
            grid->updateSources();
            solver->solve();
        }

        GLuint Eulerian2dComponent::getRenderedTexture()
        {
            renderer->draw(*grid);
            return renderer->getTextureID();
        }
    }
}
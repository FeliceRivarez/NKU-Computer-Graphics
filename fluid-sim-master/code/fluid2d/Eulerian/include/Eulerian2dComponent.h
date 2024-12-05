#pragma once
#ifndef __EULERIAN_2D_COMPONENT_H__
#define __EULERIAN_2D_COMPONENT_H__

#include "Renderer.h"
#include "Solver.h"
#include "MACGrid2d.h"

#include "Component.h"
#include "Configure.h"

#include "Logger.h"

namespace FluidSimulation {
    namespace Eulerian2d {
        class Eulerian2dComponent : public Glb::Component {
        public:
            Renderer* renderer;
            Solver* solver;
            MACGrid2d* grid;

            Eulerian2dComponent(char* description, int id) {
                this->description = description;
                this->id = id;
                renderer = NULL;
                solver = NULL;
                grid = NULL;
            }
            virtual void shutDown();
            virtual void init();
            virtual void simulate();
            virtual GLuint getRenderedTexture();
        };
    }
}


#endif
#pragma once
#ifndef __EULERIAN_3D_COMPONENT_H__
#define __EULERIAN_3D_COMPONENT_H__

#include "Renderer.h"
#include "Solver.h"
#include "MACGrid3d.h"

#include "Component.h"
#include "Configure.h"
#include "Global.h"
#include "Logger.h"

namespace FluidSimulation {
    namespace Eulerian3d {
        class Eulerian3dComponent : public Glb::Component {
        public:
            Renderer* renderer;
            Solver* solver;
            MACGrid3d* grid;

            Eulerian3dComponent(char* description, int id) {
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
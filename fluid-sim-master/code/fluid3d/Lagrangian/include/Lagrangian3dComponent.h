#pragma once
#ifndef __LAGRANGIAN_3D_COMPONENT_H__
#define __LAGRANGIAN_3D_COMPONENT_H__

#include "Renderer.h"
#include "Solver.h"
#include "ParticleSystem3d.h"

#include "Component.h"
#include "Configure.h"
#include "Logger.h"

namespace FluidSimulation
{
    namespace Lagrangian3d
    {
        class Lagrangian3dComponent : public Glb::Component
        {
        public:
            Renderer *renderer;
            Solver *solver;
            ParticleSystem3d *ps;

            Lagrangian3dComponent(char *description, int id)
            {
                this->description = description;
                this->id = id;
                renderer = NULL;
                solver = NULL;
                ps = NULL;
            }
            virtual void shutDown();
            virtual void init();
            virtual void simulate();
            virtual GLuint getRenderedTexture();
        };
    }
}

#endif
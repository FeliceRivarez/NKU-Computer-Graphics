#pragma once
#ifndef __LAGRANGIAN_3D_RENDERER_H__
#define __LAGRANGIAN_3D_RENDERER_H__

#include <glm/glm.hpp>

#include <glad/glad.h>
#include "glfw3.h"
#include "Shader.h"
#include "Container.h"

#include "Camera.h"
#include "Global.h"

#include "Configure.h"
#include <Logger.h>
#include "ParticleSystem3d.h"

namespace FluidSimulation
{

    namespace Lagrangian3d
    {

        class Renderer
        {
        public:
            Renderer(){};

            void init();

            GLuint getRenderedTexture();

            void draw(ParticleSystem3d &ps);

        private:
            void MakeVertexArrays();
        private:

            // shaders
            Glb::Shader *shader = nullptr;

            Glb::Container *container = nullptr;

            GLuint FBO = 0;
            GLuint RBO = 0;
            GLuint VAO = 0;
            GLuint VBO = 0;
            GLuint textureID = 0;

            int32_t particleNum = 0;
        };
    }
}

#endif

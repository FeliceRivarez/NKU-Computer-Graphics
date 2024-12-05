#pragma once
#ifndef __LAGRANGIAN_2D_RENDERER_H__
#define __LAGRANGIAN_2D_RENDERER_H__

#include <glad/glad.h>
#include <glfw3.h>

#include <chrono>
#include <vector>
#include "Shader.h"
#include <glm/glm.hpp>

#include "ParticleSystem2d.h"
#include "Configure.h"
#include <Logger.h>

namespace FluidSimulation
{

    namespace Lagrangian2d
    {
        class Renderer
        {
        public:
            Renderer();

            int32_t init();

            void draw(ParticleSystem2d& ps);

            GLuint getRenderedTexture();

        private:

            Glb::Shader *shader = nullptr;

            GLuint VAO = 0;
            GLuint positionVBO = 0;
            GLuint densityVBO = 0;

            GLuint FBO = 0;
            GLuint textureID = 0;
            GLuint RBO = 0;

            size_t particleNum = 0;
        };
    }
}

#endif

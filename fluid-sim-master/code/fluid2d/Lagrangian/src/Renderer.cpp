#include "Lagrangian/include/Renderer.h"

#include <iostream>
#include <fstream>
#include "Configure.h"

// This is the main definitions of renderer, which involve the OpenGL functions.
// Before browsing this file, make sure you have elementary knowledge of OpenGL.

namespace FluidSimulation
{

    namespace Lagrangian2d
    {

        Renderer::Renderer()
        {
        }

        int32_t Renderer::init()
        {
            extern std::string shaderPath;

            std::string particleVertShaderPath = shaderPath + "/DrawParticles2d.vert";
            std::string particleFragShaderPath = shaderPath + "/DrawParticles2d.frag";
            shader = new Glb::Shader();
            shader->buildFromFile(particleVertShaderPath, particleFragShaderPath);

            // generate vertex array object
            glGenVertexArrays(1, &VAO);
            // generate vertex buffer object (for position)
            glGenBuffers(1, &positionVBO);
            // generate vertex buffer object (for density)
            glGenBuffers(1, &densityVBO);

            // generate frame buffer object
            glGenFramebuffers(1, &FBO);
            // make it active
            // start fbo
            glBindFramebuffer(GL_FRAMEBUFFER, FBO);

            // generate textures
            glGenTextures(1, &textureID);
            glBindTexture(GL_TEXTURE_2D, textureID);
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, imageWidth, imageHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
            glBindTexture(GL_TEXTURE_2D, 0);

            glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, textureID, 0);

            // generate render buffer object (RBO)
            glGenRenderbuffers(1, &RBO);
            glBindRenderbuffer(GL_RENDERBUFFER, RBO);
            glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, imageWidth, imageHeight);
            glBindRenderbuffer(GL_RENDERBUFFER, 0);

            glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, RBO);
            if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
            {
                Glb::Logger::getInstance().addLog("Error: Framebuffer is not complete!");
            }
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
            // end fbo

            glViewport(0, 0, imageWidth, imageHeight);

            return 0;
        }

        void Renderer::draw(ParticleSystem2d& ps)
        {

            // bind VAO (decide which VAO we want to set)
            glBindVertexArray(VAO);

            // bind VBO to GL_ARRAY_BUFFER
            glBindBuffer(GL_ARRAY_BUFFER, positionVBO);
            // copy data to the current bound buffer(VBO)
            glBufferData(GL_ARRAY_BUFFER, ps.particles.size() * sizeof(ParticleInfo2d), ps.particles.data(), GL_STATIC_DRAW);

            glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(ParticleInfo2d), (void*)offsetof(ParticleInfo2d, position));
            // activate
            glEnableVertexAttribArray(0);

            glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(ParticleInfo2d), (void*)offsetof(ParticleInfo2d, density));
            glEnableVertexAttribArray(1);

            glBindVertexArray(0);

            particleNum = ps.particles.size();

            glBindFramebuffer(GL_FRAMEBUFFER, FBO);

            glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
            glEnable(GL_DEPTH_TEST);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            glBindVertexArray(VAO); // VAO
            shader->use();           // shader
            shader->setFloat("scale", ps.scale);

            glEnable(GL_PROGRAM_POINT_SIZE);

            glDrawArrays(GL_POINTS, 0, particleNum);

            glBindFramebuffer(GL_FRAMEBUFFER, 0);
        }

        GLuint Renderer::getRenderedTexture()
        {
            return textureID;
        }
    }

}
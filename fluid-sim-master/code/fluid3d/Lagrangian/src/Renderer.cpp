#include "fluid3d/Lagrangian/include/Renderer.h"

namespace FluidSimulation
{
    namespace Lagrangian3d
    {

        void Renderer::init()
        {
            container = new Glb::Container();
            container->resetSize(1, 1, 1);
            container->init();

            // Build Shaders

            shader = new Glb::Shader();
            std::string drawColorVertPath = shaderPath + "/DrawParticles3d.vert";
            std::string drawColorFragPath = shaderPath + "/DrawParticles3d.frag";
            shader->buildFromFile(drawColorVertPath, drawColorFragPath);

            // Generate Frame Buffers
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
            glGenBuffers(1, &VBO);
            MakeVertexArrays();

            glEnable(GL_MULTISAMPLE);

            glViewport(0, 0, imageWidth, imageHeight);
        }

        void Renderer::MakeVertexArrays()
        {
            glGenVertexArrays(1, &VAO);
            glBindVertexArray(VAO);
            glBindBuffer(GL_ARRAY_BUFFER, VBO);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(particle3d), (void *)offsetof(particle3d, position));
            glEnableVertexAttribArray(0); // location = 0
            glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(particle3d), (void *)offsetof(particle3d, density));
            glEnableVertexAttribArray(1); // location = 1
            glBindVertexArray(0);
        }

        void Renderer::draw(ParticleSystem3d &ps)
        {
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, VBO);
            glBufferData(GL_SHADER_STORAGE_BUFFER, ps.particles.size() * sizeof(particle3d), ps.particles.data(), GL_DYNAMIC_COPY);
            particleNum = ps.particles.size();

            glBindFramebuffer(GL_FRAMEBUFFER, FBO);
            glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
            glEnable(GL_DEPTH_TEST);
            glDepthFunc(GL_LEQUAL);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glEnable(GL_PROGRAM_POINT_SIZE);

            shader->use();
            shader->setMat4("view", Glb::Camera::getInstance().GetView());
            shader->setMat4("projection", Glb::Camera::getInstance().GetProjection());
            shader->setFloat("scale", ps.scale);

            glBindVertexArray(VAO);
            glDrawArrays(GL_POINTS, 0, particleNum);
            shader->unUse();

            container->draw();

            glBindFramebuffer(GL_FRAMEBUFFER, 0);
        }

        GLuint Renderer::getRenderedTexture()
        {
            return textureID;
        }
    }
}
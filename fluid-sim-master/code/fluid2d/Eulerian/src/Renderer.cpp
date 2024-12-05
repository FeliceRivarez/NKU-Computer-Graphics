#include "fluid2d/Eulerian/include/Renderer.h"
#include "stb_image.h"

namespace FluidSimulation
{
	namespace Eulerian2d
	{

		float vertices[]{
			// position	//texcood   //color
			0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f, 0.0f};

		unsigned int indices[] = {
			0, 1, 2, // first triangle
			0, 2, 3	 // second triangle
		};

		Renderer::Renderer()
		{

			extern std::string shaderPath;
			std::string vertShaderPath = shaderPath + "/DrawSmokeTexture2d.vert";
			std::string fragShaderPath = shaderPath + "/DrawSmokeTexture2d.frag";
			shader = new Glb::Shader();
			shader->buildFromFile(vertShaderPath, fragShaderPath);

			glGenVertexArrays(1, &VAO);
			glGenBuffers(1, &VBO);
			glGenBuffers(1, &EBO);

			glGenFramebuffers(1, &FBO);
			glBindFramebuffer(GL_FRAMEBUFFER, FBO);

			glGenTextures(1, &textureGridID);
			glBindTexture(GL_TEXTURE_2D, textureGridID);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, imageWidth, imageHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

			glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, textureGridID, 0);

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

			glDisable(GL_DEPTH_TEST);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

			glViewport(0, 0, imageWidth, imageHeight);

			loadTexture();
		}

		void Renderer::loadTexture()
		{
			glGenTextures(1, &smokeTexture);
			glBindTexture(GL_TEXTURE_2D, smokeTexture); // all upcoming GL_TEXTURE_2D operations now have effect on this texture object
			// set the texture wrapping parameters
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // set texture wrapping to GL_REPEAT (default wrapping method)
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
			// set texture filtering parameters
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			// load image, create texture and generate mipmaps
			int width, height, nrChannels;

			unsigned char *data = stbi_load((picturePath + "/white.png").c_str(), &width, &height, &nrChannels, 0);
			if (data)
			{
				glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
				// Mipmaps
				glGenerateMipmap(GL_TEXTURE_2D);
			}
			else
			{
				Glb::Logger::getInstance().addLog("Error: Failed to load texture!");
			}
			stbi_image_free(data);
		}

		void Renderer::draw(MACGrid2d &mGrid)
		{

			if (Eulerian2dPara::drawModel == 0)
			{
				std::vector<float> imageData;

				for (int j = 1; j <= imageHeight; j++)
				{
					for (int i = 1; i <= imageWidth; i++)
					{
						float pt_x = i * mGrid.mD.mMax[0] / (imageWidth);
						float pt_y = j * mGrid.mD.mMax[1] / (imageHeight);
						glm::vec2 pt(pt_x, pt_y);
						if (mGrid.inSolid(pt)) {
							imageData.push_back(0);
							imageData.push_back(1);
							imageData.push_back(0);
						}
						else {
							glm::vec4 color = mGrid.getRenderColor(pt);
							imageData.push_back(color.x * Eulerian2dPara::contrast);
							imageData.push_back(color.y * Eulerian2dPara::contrast);
							imageData.push_back(color.z * Eulerian2dPara::contrast);
						}
					}
				}

				glGenTextures(1, &texturePixelID);
				glBindTexture(GL_TEXTURE_2D, texturePixelID);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

				glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, imageWidth, imageHeight, 0, GL_RGB, GL_FLOAT, imageData.data());
				glBindTexture(GL_TEXTURE_2D, 0);
			}
			else
			{
				glBindFramebuffer(GL_FRAMEBUFFER, FBO);

				glClearColor(0.05f, 0.05f, 0.05f, 1.0f);
				glClear(GL_COLOR_BUFFER_BIT);

				float dt_x = mGrid.mU.mMax[0] / (Eulerian2dPara::gridNum);
				float dt_y = mGrid.mV.mMax[1] / (Eulerian2dPara::gridNum);

				for (int j = Eulerian2dPara::gridNum; j >= 1; j--)
				{
					for (int i = Eulerian2dPara::gridNum; i >= 1; i--)
					{
						float pt_x = i * mGrid.mD.mMax[0] / (Eulerian2dPara::gridNum);
						float pt_y = j * mGrid.mD.mMax[1] / (Eulerian2dPara::gridNum);

						vertices[0] = pt_x - dt_x / 2;
						vertices[1] = pt_y - dt_y / 2;
						vertices[4] = mGrid.getDensity(glm::vec2(vertices[0], vertices[1]));

						vertices[5] = pt_x + dt_x / 2;
						vertices[6] = pt_y - dt_y / 2;
						vertices[9] = mGrid.getDensity(glm::vec2(vertices[5], vertices[6]));

						vertices[10] = pt_x + dt_x / 2;
						vertices[11] = pt_y + dt_y / 2;
						vertices[14] = mGrid.getDensity(glm::vec2(vertices[10], vertices[11]));

						vertices[15] = pt_x - dt_x / 2;
						vertices[16] = pt_y + dt_y / 2;
						vertices[19] = mGrid.getDensity(glm::vec2(vertices[15], vertices[16]));

						for (int k = 0; k <= 15; k += 5)
						{
							vertices[k] = (vertices[k] / mGrid.mU.mMax[0]) * 2 - 1;
							vertices[k + 1] = (vertices[k + 1] / mGrid.mV.mMax[1]) * 2 - 1;
						}

						glBindVertexArray(VAO);
						glBindBuffer(GL_ARRAY_BUFFER, VBO);
						glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
						glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
						glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
						glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void *)0);
						glEnableVertexAttribArray(0);
						glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void *)(2 * sizeof(float)));
						glEnableVertexAttribArray(1);
						glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void *)(4 * sizeof(float)));
						glEnableVertexAttribArray(2);
						glBindVertexArray(0);

						shader->use();
						glBindTexture(GL_TEXTURE_2D, smokeTexture);
						glUniform1i(glGetUniformLocation(shader->getId(), "mTexture"), 0);
						shader->setFloat("contrast", Eulerian2dPara::contrast);
						glBindVertexArray(VAO);
						glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
					}
				}

				glBindFramebuffer(GL_FRAMEBUFFER, 0);
			}
		}

		GLuint Renderer::getTextureID()
		{
			if (Eulerian2dPara::drawModel == 0)
			{
				return texturePixelID;
			}
			return textureGridID;
		}
	}
}
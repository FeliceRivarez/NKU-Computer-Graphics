#include "fluid3d/Eulerian/include/Renderer.h"

#include "stb_image.h"

namespace FluidSimulation
{
	namespace Eulerian3d
	{

		float vertices[]{
			// position		  //texcood   //color
			0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f};

		unsigned int indices[] = {
			0, 1, 2, // first triangle
			0, 2, 3	 // second triangle
		};

		float verticesXY[] = {
			0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
			1.0f, 0.0f, 0.0f, 1.0f, 0.0f,
			1.0f, 1.0f, 0.0f, 1.0f, 1.0f,
			1.0f, 1.0f, 0.0f, 1.0f, 1.0f,
			0.0f, 1.0f, 0.0f, 0.0f, 1.0f,
			0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

		float verticesYZ[] = {
			0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
			0.0f, 1.0f, 0.0f, 1.0f, 0.0f,
			0.0f, 1.0f, 1.0f, 1.0f, 1.0f,
			0.0f, 1.0f, 1.0f, 1.0f, 1.0f,
			0.0f, 0.0f, 1.0f, 0.0f, 1.0f,
			0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

		float verticesXZ[] = {
			0.0f, 0.0f, 0.0f, 1.0f, 0.0f,
			1.0f, 0.0f, 0.0f, 0.0f, 0.0f,
			1.0f, 0.0f, 1.0f, 0.0f, 1.0f,
			1.0f, 0.0f, 1.0f, 0.0f, 1.0f,
			0.0f, 0.0f, 1.0f, 1.0f, 1.0f,
			0.0f, 0.0f, 0.0f, 1.0f, 0.0f};

		void Renderer::resetVertices(float x, float y, float z)
		{
			for (int i = 0; i < 30; i += 5) {
				if (verticesXY[i] != 0)verticesXY[i] = x;
				if (verticesXY[i+1] != 0)verticesXY[i+1] = y;
				if (verticesXY[i+2] != 0)verticesXY[i+2] = z;
				if (verticesYZ[i] != 0)verticesYZ[i] = x;
				if (verticesYZ[i + 1] != 0)verticesYZ[i + 1] = y;
				if (verticesYZ[i + 2] != 0)verticesYZ[i + 2] = z;
				if (verticesXZ[i] != 0)verticesXZ[i] = x;
				if (verticesXZ[i + 1] != 0)verticesXZ[i + 1] = y;
				if (verticesXZ[i + 2] != 0)verticesXZ[i + 2] = z;
			}
		}

		Renderer::Renderer(MACGrid3d &grid) : mGrid(grid)
		{

			container = new Glb::Container();
			float x = (float)Eulerian3dPara::theDim3d[0] / Eulerian3dPara::theDim3d[2];
			float y = (float)Eulerian3dPara::theDim3d[1] / Eulerian3dPara::theDim3d[2];
			container->resetSize(x, y, 1.0f);
			container->init();
			resetVertices(x, y, 1.0f);

			std::string particleVertShaderPath = shaderPath + "/DrawSmokePixel3d.vert";
			std::string particleFragShaderPath = shaderPath + "/DrawSmokePixel3d.frag";
			pixelShader = new Glb::Shader();
			pixelShader->buildFromFile(particleVertShaderPath, particleFragShaderPath);

			particleVertShaderPath = shaderPath + "/DrawSmokeTexture3d.vert";
			particleFragShaderPath = shaderPath + "/DrawSmokeTexture3d.frag";
			gridShader = new Glb::Shader();
			gridShader->buildFromFile(particleVertShaderPath, particleFragShaderPath);

			glGenVertexArrays(1, &VAO);
			glGenBuffers(1, &VBO);
			glGenBuffers(1, &EBO);

			glGenVertexArrays(1, &VAO_XY);
			glGenBuffers(1, &VBO_XY);
			glBindVertexArray(VAO_XY);
			glBindBuffer(GL_ARRAY_BUFFER, VBO_XY);
			glBufferData(GL_ARRAY_BUFFER, sizeof(verticesXY), verticesXY, GL_STATIC_DRAW);
			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void *)0);
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void *)(3 * sizeof(float)));
			glEnableVertexAttribArray(1);
			glBindVertexArray(0);

			glGenVertexArrays(1, &VAO_YZ);
			glGenBuffers(1, &VBO_YZ);
			glBindVertexArray(VAO_YZ);
			glBindBuffer(GL_ARRAY_BUFFER, VBO_YZ);
			glBufferData(GL_ARRAY_BUFFER, sizeof(verticesYZ), verticesYZ, GL_STATIC_DRAW);
			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void *)0);
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void *)(3 * sizeof(float)));
			glEnableVertexAttribArray(1);
			glBindVertexArray(0);

			glGenVertexArrays(1, &VAO_XZ);
			glGenBuffers(1, &VBO_XZ);
			glBindVertexArray(VAO_XZ);
			glBindBuffer(GL_ARRAY_BUFFER, VBO_XZ);
			glBufferData(GL_ARRAY_BUFFER, sizeof(verticesXZ), verticesXZ, GL_STATIC_DRAW);
			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void *)0);
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void *)(3 * sizeof(float)));
			glEnableVertexAttribArray(1);
			glBindVertexArray(0);

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

			// FBO
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

			glDisable(GL_DEPTH_TEST);

			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

			dataXY = new float[4 * pixelX * pixelY];
			dataYZ = new float[4 * pixelY * pixelZ];
			dataXZ = new float[4 * pixelX * pixelZ];

			loadTexture();

			glViewport(0, 0, imageWidth, imageHeight);
		}

		void Renderer::draw()
		{
			glBindFramebuffer(GL_FRAMEBUFFER, FBO);

			glClearColor(0.05f, 0.05f, 0.05f, 1.0f);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			if (Eulerian3dPara::oneSheet)
			{
				drawOneSheet();
			}
			else
			{
				if (Eulerian3dPara::xySheetsON)
				{
					drawXYSheets();
				}
				if (Eulerian3dPara::yzSheetsON)
				{
					drawYZSheets();
				}
				if (Eulerian3dPara::xzSheetsON)
				{
					drawXZSheets();
				}
			}
			container->draw();
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
		}

		void Renderer::drawOneSheet()
		{
			if (Eulerian3dPara::xySheetsON)
			{
				drawOneSheetXY();
			}
			if (Eulerian3dPara::yzSheetsON)
			{
				drawOneSheetYZ();
			}
			if (Eulerian3dPara::xzSheetsON)
			{
				drawOneSheetXZ();
			}

		}

		void Renderer::drawOneSheetXY()
		{
			if (Eulerian3dPara::drawModel == 0)
			{
				for (int i = 1; i <= pixelX; i++)
				{
					for (int j = 1; j <= pixelY; j++)
					{
						float pt_x = i * mGrid.mD.mMax[0] / (pixelX);
						float pt_y = j * mGrid.mD.mMax[1] / (pixelY);
						float pt_z = Eulerian3dPara::distanceZ * mGrid.mD.mMax[2];
						glm::vec3 pt(pt_x, pt_y, pt_z);
						glm::vec4 color = mGrid.getRenderColor(pt);
						dataXY[4 * ((pixelY - j) * pixelX + (i - 1))] = color.r;
						dataXY[4 * ((pixelY - j) * pixelX + (i - 1)) + 1] = color.g;
						dataXY[4 * ((pixelY - j) * pixelX + (i - 1)) + 2] = color.b;
						dataXY[4 * ((pixelY - j) * pixelX + (i - 1)) + 3] = color.a;
					}
				}

				unsigned int texture;
				glGenTextures(1, &texture);
				glBindTexture(GL_TEXTURE_2D, texture);
				// set the texture wrapping parameters
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
				// set texture filtering parameters
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
				glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, pixelX, pixelY, 0, GL_RGBA, GL_FLOAT, dataXY);
				glGenerateMipmap(GL_TEXTURE_2D);

				pixelShader->use();
				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, texture);
				glUniform1i(glGetUniformLocation(pixelShader->getId(), "aTexture"), 0);

				glm::mat4 view = Glb::Camera::getInstance().GetView();
				glm::mat4 projection = Glb::Camera::getInstance().GetProjection();

				glm::mat4 model = glm::mat4(1.0f);
				model = glm::translate(model, glm::vec3(0.0f, 0.0f, Eulerian3dPara::distanceZ));
				glUniformMatrix4fv(glGetUniformLocation(pixelShader->getId(), "view"), 1, GL_FALSE, glm::value_ptr(view));
				glUniformMatrix4fv(glGetUniformLocation(pixelShader->getId(), "projection"), 1, GL_FALSE, glm::value_ptr(projection));
				glUniformMatrix4fv(glGetUniformLocation(pixelShader->getId(), "model"), 1, GL_FALSE, glm::value_ptr(model));
				pixelShader->setFloat("contrast", Eulerian3dPara::contrast);
				glBindVertexArray(VAO_XY);
				pixelShader->use();
				glDrawArrays(GL_TRIANGLES, 0, 6);
			}
			else if (Eulerian3dPara::drawModel == 1)
			{
				float dt_x = mGrid.mD.mMax[0] / (Eulerian3dPara::gridNumX);
				float dt_y = mGrid.mD.mMax[1] / (Eulerian3dPara::gridNumY);

				for (int i = 1; i <= Eulerian3dPara::gridNumX; i++)
				{
					for (int j = 1; j <= Eulerian3dPara::gridNumY; j++)
					{
						float pt_x = i * mGrid.mD.mMax[0] / (Eulerian3dPara::gridNumX);
						float pt_y = j * mGrid.mD.mMax[1] / (Eulerian3dPara::gridNumY);
						float pt_z = Eulerian3dPara::distanceZ * mGrid.mD.mMax[2];

						vertices[0] = pt_x + dt_x / 2;
						vertices[1] = pt_y - dt_y / 2;
						vertices[5] = mGrid.getDensity(glm::vec3(vertices[0], vertices[1], pt_z));

						vertices[6] = pt_x + dt_x / 2;
						vertices[7] = pt_y + dt_y / 2;
						vertices[11] = mGrid.getDensity(glm::vec3(vertices[6], vertices[7], pt_z));

						vertices[12] = pt_x - dt_x / 2;
						vertices[13] = pt_y + dt_y / 2;
						vertices[17] = mGrid.getDensity(glm::vec3(vertices[12], vertices[13], pt_z));

						vertices[18] = pt_x - dt_x / 2;
						vertices[19] = pt_y - dt_y / 2;
						vertices[23] = mGrid.getDensity(glm::vec3(vertices[18], vertices[19], pt_z));

						for (int k = 0; k <= 18; k += 6)
						{
							vertices[k] = (vertices[k] / mGrid.mD.mMax[2]);
							vertices[k + 1] = (vertices[k + 1] / mGrid.mD.mMax[2]);
							vertices[k + 2] = 0;
						}

						glBindVertexArray(VAO);
						glBindBuffer(GL_ARRAY_BUFFER, VBO);
						glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
						glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
						glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
						glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)0);
						glEnableVertexAttribArray(0);
						glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)(3 * sizeof(float)));
						glEnableVertexAttribArray(1);
						glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)(5 * sizeof(float)));
						glEnableVertexAttribArray(2);
						glBindVertexArray(0);

						gridShader->use();
						glBindTexture(GL_TEXTURE_2D, smokeTexture);
						glUniform1i(glGetUniformLocation(gridShader->getId(), "mTexture"), 0);

						glm::mat4 view = Glb::Camera::getInstance().GetView();
						glm::mat4 projection = Glb::Camera::getInstance().GetProjection();

						glm::mat4 model = glm::mat4(1.0f);
						model = glm::translate(model, glm::vec3(0.0f, 0.0f, Eulerian3dPara::distanceZ));
						glUniformMatrix4fv(glGetUniformLocation(gridShader->getId(), "view"), 1, GL_FALSE, glm::value_ptr(view));
						glUniformMatrix4fv(glGetUniformLocation(gridShader->getId(), "projection"), 1, GL_FALSE, glm::value_ptr(projection));
						glUniformMatrix4fv(glGetUniformLocation(gridShader->getId(), "model"), 1, GL_FALSE, glm::value_ptr(model));
						gridShader->setFloat("contrast", Eulerian3dPara::contrast);

						glBindVertexArray(VAO);
						glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
					}
				}
			}
		}

		void Renderer::drawOneSheetXZ()
		{
			if (Eulerian3dPara::drawModel == 0)
			{
				for (int k = pixelZ; k >= 1; k--)
				{
					for (int i = pixelX; i >= 1; i--)
					{
						float pt_x = i * mGrid.mD.mMax[0] / (pixelX);
						float pt_y = Eulerian3dPara::distanceY * mGrid.mD.mMax[1];
						float pt_z = k * mGrid.mD.mMax[2] / (pixelZ);
						glm::vec3 pt(pt_x, pt_y, pt_z);
						glm::vec4 color = mGrid.getRenderColor(pt);
						dataXZ[4 * ((pixelZ - k) * pixelX + (pixelX - i))] = color.r;
						dataXZ[4 * ((pixelZ - k) * pixelX + (pixelX - i)) + 1] = color.g;
						dataXZ[4 * ((pixelZ - k) * pixelX + (pixelX - i)) + 2] = color.b;
						dataXZ[4 * ((pixelZ - k) * pixelX + (pixelX - i)) + 3] = color.a;
					}
				}

				unsigned int texture;
				glGenTextures(1, &texture);
				glBindTexture(GL_TEXTURE_2D, texture); // all upcoming GL_TEXTURE_2D operations now have effect on this texture object
				// set the texture wrapping parameters
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // set texture wrapping to GL_REPEAT (default wrapping method)
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
				// set texture filtering parameters
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
				glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, pixelX, pixelZ, 0, GL_RGBA, GL_FLOAT, dataXZ);
				glGenerateMipmap(GL_TEXTURE_2D);

				pixelShader->use();
				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, texture);
				glUniform1i(glGetUniformLocation(pixelShader->getId(), "aTexture"), 0);

				glm::mat4 view = Glb::Camera::getInstance().GetView();
				glm::mat4 projection = Glb::Camera::getInstance().GetProjection();

				glm::mat4 model = glm::mat4(1.0f);
				model = glm::translate(model, glm::vec3(0.0f, Eulerian3dPara::distanceY, 0.0f));
				glUniformMatrix4fv(glGetUniformLocation(pixelShader->getId(), "view"), 1, GL_FALSE, glm::value_ptr(view));
				glUniformMatrix4fv(glGetUniformLocation(pixelShader->getId(), "projection"), 1, GL_FALSE, glm::value_ptr(projection));
				glUniformMatrix4fv(glGetUniformLocation(pixelShader->getId(), "model"), 1, GL_FALSE, glm::value_ptr(model));
				pixelShader->setFloat("contrast", Eulerian3dPara::contrast);
				glBindVertexArray(VAO_XZ);
				pixelShader->use();
				glDrawArrays(GL_TRIANGLES, 0, 6);
			}
			else if (Eulerian3dPara::drawModel == 1)
			{

				float dt_x = mGrid.mD.mMax[0] / (Eulerian3dPara::gridNumX);
				float dt_z = mGrid.mD.mMax[2] / (Eulerian3dPara::gridNumZ);

				for (int k = Eulerian3dPara::gridNumZ; k >= 1; k--)
				{
					for (int i = Eulerian3dPara::gridNumX; i >= 1; i--)
					{
						float pt_x = i * mGrid.mD.mMax[0] / (Eulerian3dPara::gridNumX);
						float pt_y = Eulerian3dPara::distanceY * mGrid.mD.mMax[1];
						float pt_z = k * mGrid.mD.mMax[2] / (Eulerian3dPara::gridNumZ);

						vertices[0] = pt_x + dt_x / 2;
						vertices[2] = pt_z - dt_z / 2;
						vertices[5] = mGrid.getDensity(glm::vec3(vertices[0], pt_y, vertices[2]));

						vertices[6] = pt_x - dt_x / 2;
						vertices[8] = pt_z - dt_z / 2;
						vertices[11] = mGrid.getDensity(glm::vec3(vertices[6], pt_y, vertices[8]));

						vertices[12] = pt_x - dt_x / 2;
						vertices[14] = pt_z + dt_z / 2;
						vertices[17] = mGrid.getDensity(glm::vec3(vertices[12], pt_y, vertices[14]));

						vertices[18] = pt_x + dt_x / 2;
						vertices[20] = pt_z + dt_z / 2;
						vertices[23] = mGrid.getDensity(glm::vec3(vertices[18], pt_y, vertices[20]));

						for (int k = 0; k <= 18; k += 6)
						{
							vertices[k] = (vertices[k] / mGrid.mD.mMax[2]);
							vertices[k + 1] = 0;
							vertices[k + 2] = (vertices[k + 2] / mGrid.mD.mMax[2]);
						}

						glBindVertexArray(VAO);
						glBindBuffer(GL_ARRAY_BUFFER, VBO);
						glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
						glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
						glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
						glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)0);
						glEnableVertexAttribArray(0);
						glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)(3 * sizeof(float)));
						glEnableVertexAttribArray(1);
						glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)(5 * sizeof(float)));
						glEnableVertexAttribArray(2);
						glBindVertexArray(0);

						gridShader->use();
						glBindTexture(GL_TEXTURE_2D, smokeTexture);
						glUniform1i(glGetUniformLocation(gridShader->getId(), "mTexture"), 0);

						glm::mat4 view = Glb::Camera::getInstance().GetView();
						glm::mat4 projection = Glb::Camera::getInstance().GetProjection();

						glm::mat4 model = glm::mat4(1.0f);
						model = glm::translate(model, glm::vec3(0.0f, Eulerian3dPara::distanceY * mGrid.mD.mMax[1] / mGrid.mD.mMax[2], 0.0f));
						glUniformMatrix4fv(glGetUniformLocation(gridShader->getId(), "view"), 1, GL_FALSE, glm::value_ptr(view));
						glUniformMatrix4fv(glGetUniformLocation(gridShader->getId(), "projection"), 1, GL_FALSE, glm::value_ptr(projection));
						glUniformMatrix4fv(glGetUniformLocation(gridShader->getId(), "model"), 1, GL_FALSE, glm::value_ptr(model));
						gridShader->setFloat("contrast", Eulerian3dPara::contrast);
						glBindVertexArray(VAO);
						glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
					}
				}
			}
		}

		void Renderer::drawOneSheetYZ()
		{
			if (Eulerian3dPara::drawModel == 0)
			{

				for (int k = pixelZ; k >= 1; k--)
				{
					for (int j = 1; j <= pixelY; j++)
					{
						float pt_x = Eulerian3dPara::distanceX * mGrid.mD.mMax[0];
						float pt_y = j * mGrid.mD.mMax[1] / (pixelY);
						float pt_z = k * mGrid.mD.mMax[2] / (pixelZ);
						glm::vec3 pt(pt_x, pt_y, pt_z);
						glm::vec4 color = mGrid.getRenderColor(pt);
						dataYZ[4 * ((pixelZ - k) * pixelY + (j - 1))] = color.r;
						dataYZ[4 * ((pixelZ - k) * pixelY + (j - 1)) + 1] = color.g;
						dataYZ[4 * ((pixelZ - k) * pixelY + (j - 1)) + 2] = color.b;
						dataYZ[4 * ((pixelZ - k) * pixelY + (j - 1)) + 3] = color.a;
					}
				}

				unsigned int texture;
				glGenTextures(1, &texture);
				glBindTexture(GL_TEXTURE_2D, texture); // all upcoming GL_TEXTURE_2D operations now have effect on this texture object
				// set the texture wrapping parameters
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // set texture wrapping to GL_REPEAT (default wrapping method)
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
				// set texture filtering parameters
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
				glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
				glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, pixelY, pixelZ, 0, GL_RGBA, GL_FLOAT, dataYZ);
				glGenerateMipmap(GL_TEXTURE_2D);

				pixelShader->use();
				glActiveTexture(GL_TEXTURE0);
				glBindTexture(GL_TEXTURE_2D, texture);
				glUniform1i(glGetUniformLocation(pixelShader->getId(), "aTexture"), 0);

				glm::mat4 view = Glb::Camera::getInstance().GetView();
				glm::mat4 projection = Glb::Camera::getInstance().GetProjection();

				glm::mat4 model = glm::mat4(1.0f);
				model = glm::translate(model, glm::vec3((Eulerian3dPara::distanceX * mGrid.mD.mMax[0] / mGrid.mD.mMax[2]), 0.0f, 0.0f));
				glUniformMatrix4fv(glGetUniformLocation(pixelShader->getId(), "view"), 1, GL_FALSE, glm::value_ptr(view));
				glUniformMatrix4fv(glGetUniformLocation(pixelShader->getId(), "projection"), 1, GL_FALSE, glm::value_ptr(projection));
				glUniformMatrix4fv(glGetUniformLocation(pixelShader->getId(), "model"), 1, GL_FALSE, glm::value_ptr(model));
				pixelShader->setFloat("contrast", Eulerian3dPara::contrast);
				glBindVertexArray(VAO_YZ);
				pixelShader->use();
				glDrawArrays(GL_TRIANGLES, 0, 6);
			}
			else if (Eulerian3dPara::drawModel == 1)
			{

				float dt_y = mGrid.mD.mMax[1] / (Eulerian3dPara::gridNumY);
				float dt_z = mGrid.mD.mMax[2] / (Eulerian3dPara::gridNumZ);

				for (int k = Eulerian3dPara::gridNumZ; k >= 1; k--)
				{
					for (int j = 1; j <= Eulerian3dPara::gridNumY; j++)
					{
						float pt_x = Eulerian3dPara::distanceX * mGrid.mD.mMax[0];
						float pt_y = j * mGrid.mD.mMax[1] / (Eulerian3dPara::gridNumY);
						float pt_z = k * mGrid.mD.mMax[2] / (Eulerian3dPara::gridNumZ);

						vertices[1] = pt_y - dt_y / 2;
						vertices[2] = pt_z - dt_z / 2;
						vertices[5] = mGrid.getDensity(glm::vec3(pt_x, vertices[1], vertices[2]));

						vertices[7] = pt_y + dt_y / 2;
						vertices[8] = pt_z - dt_z / 2;
						vertices[11] = mGrid.getDensity(glm::vec3(pt_x, vertices[7], vertices[8]));

						vertices[13] = pt_y + dt_y / 2;
						vertices[14] = pt_z + dt_z / 2;
						vertices[17] = mGrid.getDensity(glm::vec3(pt_x, vertices[13], vertices[14]));

						vertices[19] = pt_y - dt_y / 2;
						vertices[20] = pt_z + dt_z / 2;
						vertices[23] = mGrid.getDensity(glm::vec3(pt_x, vertices[19], vertices[20]));

						for (int k = 1; k <= 19; k += 6)
						{
							vertices[k - 1] = 0;
							vertices[k] = (vertices[k] / mGrid.mD.mMax[1]);
							vertices[k + 1] = (vertices[k + 1] / mGrid.mD.mMax[2]);
						}

						glBindVertexArray(VAO);
						glBindBuffer(GL_ARRAY_BUFFER, VBO);
						glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
						glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
						glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
						glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)0);
						glEnableVertexAttribArray(0);
						glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)(3 * sizeof(float)));
						glEnableVertexAttribArray(1);
						glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)(5 * sizeof(float)));
						glEnableVertexAttribArray(2);
						glBindVertexArray(0);

						gridShader->use();
						glBindTexture(GL_TEXTURE_2D, smokeTexture);
						glUniform1i(glGetUniformLocation(gridShader->getId(), "mTexture"), 0);

						glm::mat4 view = Glb::Camera::getInstance().GetView();
						glm::mat4 projection = Glb::Camera::getInstance().GetProjection();

						glm::mat4 model = glm::mat4(1.0f);
						model = glm::translate(model, glm::vec3(Eulerian3dPara::distanceX * mGrid.mD.mMax[0] / mGrid.mD.mMax[2], 0.0f, 0.0f));
						glUniformMatrix4fv(glGetUniformLocation(gridShader->getId(), "view"), 1, GL_FALSE, glm::value_ptr(view));
						glUniformMatrix4fv(glGetUniformLocation(gridShader->getId(), "projection"), 1, GL_FALSE, glm::value_ptr(projection));
						glUniformMatrix4fv(glGetUniformLocation(gridShader->getId(), "model"), 1, GL_FALSE, glm::value_ptr(model));
						gridShader->setFloat("contrast", Eulerian3dPara::contrast);
						glBindVertexArray(VAO);
						glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
					}
				}
			}
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

		void Renderer::drawXYSheets()
		{
			if (Eulerian3dPara::drawModel == 0)
			{

				for (float k = eps; k <= mGrid.mD.mMax[2]; k += (mGrid.mD.mMax[2] - 2 * eps) / (Eulerian3dPara::xySheetsNum - 1))
				{

					for (int i = 1; i <= pixelX; i++)
					{
						for (int j = 1; j <= pixelY; j++)
						{
							float pt_x = i * mGrid.mD.mMax[0] / (pixelX);
							float pt_y = j * mGrid.mD.mMax[1] / (pixelY);
							float pt_z = k;
							glm::vec3 pt(pt_x, pt_y, pt_z);
							glm::vec4 color = mGrid.getRenderColor(pt);
							dataXY[4 * ((pixelY - j) * pixelX + (i - 1))] = color.r;
							dataXY[4 * ((pixelY - j) * pixelX + (i - 1)) + 1] = color.g;
							dataXY[4 * ((pixelY - j) * pixelX + (i - 1)) + 2] = color.b;
							dataXY[4 * ((pixelY - j) * pixelX + (i - 1)) + 3] = color.a;
						}
					}

					unsigned int texture;
					glGenTextures(1, &texture);
					glBindTexture(GL_TEXTURE_2D, texture);
					// set the texture wrapping parameters
					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
					// set texture filtering parameters
					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
					glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, pixelX, pixelY, 0, GL_RGBA, GL_FLOAT, dataXY);
					glGenerateMipmap(GL_TEXTURE_2D);

					pixelShader->use();
					glActiveTexture(GL_TEXTURE0);
					glBindTexture(GL_TEXTURE_2D, texture);
					glUniform1i(glGetUniformLocation(pixelShader->getId(), "aTexture"), 0);

					glm::mat4 view = Glb::Camera::getInstance().GetView();
					glm::mat4 projection = Glb::Camera::getInstance().GetProjection();

					glm::mat4 model = glm::mat4(1.0f);
					model = glm::translate(model, glm::vec3(0.0f, 0.0f, k / mGrid.mD.mMax[2]));
					glUniformMatrix4fv(glGetUniformLocation(pixelShader->getId(), "view"), 1, GL_FALSE, glm::value_ptr(view));
					glUniformMatrix4fv(glGetUniformLocation(pixelShader->getId(), "projection"), 1, GL_FALSE, glm::value_ptr(projection));
					glUniformMatrix4fv(glGetUniformLocation(pixelShader->getId(), "model"), 1, GL_FALSE, glm::value_ptr(model));
					pixelShader->setFloat("contrast", Eulerian3dPara::contrast);
					glBindVertexArray(VAO_XY);
					pixelShader->use();
					glDrawArrays(GL_TRIANGLES, 0, 6);
				}
			}
			else if (Eulerian3dPara::drawModel == 1)
			{
				float dt_x = mGrid.mD.mMax[0] / (Eulerian3dPara::gridNumX);
				float dt_y = mGrid.mD.mMax[1] / (Eulerian3dPara::gridNumY);
				float dt_z = 1.0f / (Eulerian3dPara::xySheetsNum + 1);

				for (float k = dt_z; k < 1; k += dt_z)
				{
					for (int i = 1; i <= Eulerian3dPara::gridNumX; i++)
					{
						for (int j = 1; j <= Eulerian3dPara::gridNumY; j++)
						{
							float pt_x = i * mGrid.mD.mMax[0] / (Eulerian3dPara::gridNumX);
							float pt_y = j * mGrid.mD.mMax[1] / (Eulerian3dPara::gridNumY);
							float pt_z = k * mGrid.mD.mMax[2];

							vertices[0] = pt_x + dt_x / 2;
							vertices[1] = pt_y - dt_y / 2;
							vertices[5] = mGrid.getDensity(glm::vec3(vertices[0], vertices[1], pt_z));

							vertices[6] = pt_x + dt_x / 2;
							vertices[7] = pt_y + dt_y / 2;
							vertices[11] = mGrid.getDensity(glm::vec3(vertices[6], vertices[7], pt_z));

							vertices[12] = pt_x - dt_x / 2;
							vertices[13] = pt_y + dt_y / 2;
							vertices[17] = mGrid.getDensity(glm::vec3(vertices[12], vertices[13], pt_z));

							vertices[18] = pt_x - dt_x / 2;
							vertices[19] = pt_y - dt_y / 2;
							vertices[23] = mGrid.getDensity(glm::vec3(vertices[18], vertices[19], pt_z));

							for (int k = 0; k <= 18; k += 6)
							{
								vertices[k] = (vertices[k] / mGrid.mD.mMax[2]);
								vertices[k + 1] = (vertices[k + 1] / mGrid.mD.mMax[2]);
								vertices[k + 2] = 0;
							}

							glBindVertexArray(VAO);
							glBindBuffer(GL_ARRAY_BUFFER, VBO);
							glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
							glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
							glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
							glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)0);
							glEnableVertexAttribArray(0);
							glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)(3 * sizeof(float)));
							glEnableVertexAttribArray(1);
							glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)(5 * sizeof(float)));
							glEnableVertexAttribArray(2);
							glBindVertexArray(0);

							gridShader->use();
							glBindTexture(GL_TEXTURE_2D, smokeTexture);
							glUniform1i(glGetUniformLocation(gridShader->getId(), "mTexture"), 0);

							glm::mat4 view = Glb::Camera::getInstance().GetView();
							glm::mat4 projection = Glb::Camera::getInstance().GetProjection();

							glm::mat4 model = glm::mat4(1.0f);
							model = glm::translate(model, glm::vec3(0.0f, 0.0f, k));
							glUniformMatrix4fv(glGetUniformLocation(gridShader->getId(), "view"), 1, GL_FALSE, glm::value_ptr(view));
							glUniformMatrix4fv(glGetUniformLocation(gridShader->getId(), "projection"), 1, GL_FALSE, glm::value_ptr(projection));
							glUniformMatrix4fv(glGetUniformLocation(gridShader->getId(), "model"), 1, GL_FALSE, glm::value_ptr(model));
							gridShader->setFloat("contrast", Eulerian3dPara::contrast);
							glBindVertexArray(VAO);
							glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
						}
					}
				}
			}
		}

		void Renderer::drawYZSheets()
		{
			if (Eulerian3dPara::drawModel == 0)
			{
				for (float i = eps; i <= mGrid.mD.mMax[0]; i += (mGrid.mD.mMax[0] - 2 * eps) / (Eulerian3dPara::yzSheetsNum - 1))
				{
					for (int k = pixelZ; k >= 1; k--)
					{
						for (int j = 1; j <= pixelY; j++)
						{
							float pt_x = i;
							float pt_y = j * mGrid.mD.mMax[1] / (pixelY);
							float pt_z = k * mGrid.mD.mMax[2] / (pixelZ);
							glm::vec3 pt(pt_x, pt_y, pt_z);
							glm::vec4 color = mGrid.getRenderColor(pt);
							dataYZ[4 * ((pixelZ - k) * pixelY + (j - 1))] = color.r;
							dataYZ[4 * ((pixelZ - k) * pixelY + (j - 1)) + 1] = color.g;
							dataYZ[4 * ((pixelZ - k) * pixelY + (j - 1)) + 2] = color.b;
							dataYZ[4 * ((pixelZ - k) * pixelY + (j - 1)) + 3] = color.a;
						}
					}

					unsigned int texture;
					glGenTextures(1, &texture);
					glBindTexture(GL_TEXTURE_2D, texture); // all upcoming GL_TEXTURE_2D operations now have effect on this texture object
					// set the texture wrapping parameters
					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // set texture wrapping to GL_REPEAT (default wrapping method)
					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
					// set texture filtering parameters
					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
					glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, pixelY, pixelZ, 0, GL_RGBA, GL_FLOAT, dataYZ);
					glGenerateMipmap(GL_TEXTURE_2D);

					pixelShader->use();
					glActiveTexture(GL_TEXTURE0);
					glBindTexture(GL_TEXTURE_2D, texture);
					glUniform1i(glGetUniformLocation(pixelShader->getId(), "aTexture"), 0);

					glm::mat4 view = Glb::Camera::getInstance().GetView();
					glm::mat4 projection = Glb::Camera::getInstance().GetProjection();

					glm::mat4 model = glm::mat4(1.0f);
					model = glm::translate(model, glm::vec3(i / mGrid.mD.mMax[2], 0.0f, 0.0f));
					glUniformMatrix4fv(glGetUniformLocation(pixelShader->getId(), "view"), 1, GL_FALSE, glm::value_ptr(view));
					glUniformMatrix4fv(glGetUniformLocation(pixelShader->getId(), "projection"), 1, GL_FALSE, glm::value_ptr(projection));
					glUniformMatrix4fv(glGetUniformLocation(pixelShader->getId(), "model"), 1, GL_FALSE, glm::value_ptr(model));
					pixelShader->setFloat("contrast", Eulerian3dPara::contrast);
					glBindVertexArray(VAO_YZ);
					pixelShader->use();
					glDrawArrays(GL_TRIANGLES, 0, 6);
				}
			}
			else if (Eulerian3dPara::drawModel == 1)
			{

				float dt_y = mGrid.mD.mMax[1] / (Eulerian3dPara::gridNumY);
				float dt_z = mGrid.mD.mMax[2] / (Eulerian3dPara::gridNumZ);

				float dt_x = 1.0f / (Eulerian3dPara::yzSheetsNum + 1);
				for (float i = dt_x; i < 0.999f; i += dt_x)
				{
					for (int k = Eulerian3dPara::gridNumZ; k >= 1; k--)
					{
						for (int j = 1; j <= Eulerian3dPara::gridNumY; j++)
						{
							float pt_x = i * mGrid.mD.mMax[0];
							float pt_y = j * mGrid.mD.mMax[1] / (Eulerian3dPara::gridNumY);
							float pt_z = k * mGrid.mD.mMax[2] / (Eulerian3dPara::gridNumZ);

							vertices[1] = pt_y - dt_y / 2;
							vertices[2] = pt_z - dt_z / 2;
							vertices[5] = mGrid.getDensity(glm::vec3(pt_x, vertices[1], vertices[2]));

							vertices[7] = pt_y + dt_y / 2;
							vertices[8] = pt_z - dt_z / 2;
							vertices[11] = mGrid.getDensity(glm::vec3(pt_x, vertices[7], vertices[8]));

							vertices[13] = pt_y + dt_y / 2;
							vertices[14] = pt_z + dt_z / 2;
							vertices[17] = mGrid.getDensity(glm::vec3(pt_x, vertices[13], vertices[14]));

							vertices[19] = pt_y - dt_y / 2;
							vertices[20] = pt_z + dt_z / 2;
							vertices[23] = mGrid.getDensity(glm::vec3(pt_x, vertices[19], vertices[20]));

							for (int k = 1; k <= 19; k += 6)
							{
								vertices[k - 1] = 0;
								vertices[k] = (vertices[k] / mGrid.mD.mMax[2]);
								vertices[k + 1] = (vertices[k + 1] / mGrid.mD.mMax[2]);
							}

							glBindVertexArray(VAO);
							glBindBuffer(GL_ARRAY_BUFFER, VBO);
							glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
							glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
							glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
							glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)0);
							glEnableVertexAttribArray(0);
							glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)(3 * sizeof(float)));
							glEnableVertexAttribArray(1);
							glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)(5 * sizeof(float)));
							glEnableVertexAttribArray(2);
							glBindVertexArray(0);

							gridShader->use();
							glBindTexture(GL_TEXTURE_2D, smokeTexture);
							glUniform1i(glGetUniformLocation(gridShader->getId(), "mTexture"), 0);

							glm::mat4 view = Glb::Camera::getInstance().GetView();
							glm::mat4 projection = Glb::Camera::getInstance().GetProjection();

							glm::mat4 model = glm::mat4(1.0f);
							model = glm::translate(model, glm::vec3(i * mGrid.mD.mMax[0] / mGrid.mD.mMax[2], 0.0f, 0.0f));
							glUniformMatrix4fv(glGetUniformLocation(gridShader->getId(), "view"), 1, GL_FALSE, glm::value_ptr(view));
							glUniformMatrix4fv(glGetUniformLocation(gridShader->getId(), "projection"), 1, GL_FALSE, glm::value_ptr(projection));
							glUniformMatrix4fv(glGetUniformLocation(gridShader->getId(), "model"), 1, GL_FALSE, glm::value_ptr(model));
							gridShader->setFloat("contrast", Eulerian3dPara::contrast);
							glBindVertexArray(VAO);
							glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
						}
					}
				}
			}
		}

		void Renderer::drawXZSheets()
		{
			if (Eulerian3dPara::drawModel == 0)
			{

				for (float j = eps; j <= mGrid.mD.mMax[1]; j += (mGrid.mD.mMax[1] - 2 * eps) / (Eulerian3dPara::xzSheetsNum - 1))
				{
					for (int k = pixelZ; k >= 1; k--)
					{
						for (int i = pixelX; i >= 1; i--)
						{
							float pt_x = i * mGrid.mD.mMax[0] / (pixelX);
							float pt_y = j;
							float pt_z = k * mGrid.mD.mMax[2] / (pixelZ);
							glm::vec3 pt(pt_x, pt_y, pt_z);
							glm::vec4 color = mGrid.getRenderColor(pt);
							dataXZ[4 * ((pixelZ - k) * pixelX + (pixelX - i))] = color.r;
							dataXZ[4 * ((pixelZ - k) * pixelX + (pixelX - i)) + 1] = color.g;
							dataXZ[4 * ((pixelZ - k) * pixelX + (pixelX - i)) + 2] = color.b;
							dataXZ[4 * ((pixelZ - k) * pixelX + (pixelX - i)) + 3] = color.a;
						}
					}

					unsigned int texture;
					glGenTextures(1, &texture);
					glBindTexture(GL_TEXTURE_2D, texture); // all upcoming GL_TEXTURE_2D operations now have effect on this texture object
					// set the texture wrapping parameters
					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // set texture wrapping to GL_REPEAT (default wrapping method)
					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
					// set texture filtering parameters
					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
					glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
					glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, pixelX, pixelZ, 0, GL_RGBA, GL_FLOAT, dataXZ);
					glGenerateMipmap(GL_TEXTURE_2D);

					pixelShader->use();
					glActiveTexture(GL_TEXTURE0);
					glBindTexture(GL_TEXTURE_2D, texture);
					glUniform1i(glGetUniformLocation(pixelShader->getId(), "aTexture"), 0);

					glm::mat4 view = Glb::Camera::getInstance().GetView();
					glm::mat4 projection = Glb::Camera::getInstance().GetProjection();

					glm::mat4 model = glm::mat4(1.0f);
					model = glm::translate(model, glm::vec3(0.0f, j / mGrid.mD.mMax[1], 0.0f));
					glUniformMatrix4fv(glGetUniformLocation(pixelShader->getId(), "view"), 1, GL_FALSE, glm::value_ptr(view));
					glUniformMatrix4fv(glGetUniformLocation(pixelShader->getId(), "projection"), 1, GL_FALSE, glm::value_ptr(projection));
					glUniformMatrix4fv(glGetUniformLocation(pixelShader->getId(), "model"), 1, GL_FALSE, glm::value_ptr(model));
					pixelShader->setFloat("contrast", Eulerian3dPara::contrast);
					glBindVertexArray(VAO_XZ);
					pixelShader->use();
					glDrawArrays(GL_TRIANGLES, 0, 6);
				}
			}
			else if (Eulerian3dPara::drawModel == 1)
			{

				float dt_x = mGrid.mD.mMax[0] / (Eulerian3dPara::gridNumX);
				float dt_z = mGrid.mD.mMax[2] / (Eulerian3dPara::gridNumZ);

				float dt_y = 1.0f / (Eulerian3dPara::xzSheetsNum + 1);

				for (float j = dt_y; j < 1; j += dt_y)
				{

					for (int k = Eulerian3dPara::gridNumZ; k >= 1; k--)
					{
						for (int i = Eulerian3dPara::gridNumX; i >= 1; i--)
						{
							float pt_x = i * mGrid.mD.mMax[0] / (Eulerian3dPara::gridNumX);
							float pt_y = j * mGrid.mD.mMax[1];
							float pt_z = k * mGrid.mD.mMax[2] / (Eulerian3dPara::gridNumZ);

							vertices[0] = pt_x + dt_x / 2;
							vertices[2] = pt_z - dt_z / 2;
							vertices[5] = mGrid.getDensity(glm::vec3(vertices[0], pt_y, vertices[2]));

							vertices[6] = pt_x - dt_x / 2;
							vertices[8] = pt_z - dt_z / 2;
							vertices[11] = mGrid.getDensity(glm::vec3(vertices[6], pt_y, vertices[8]));

							vertices[12] = pt_x - dt_x / 2;
							vertices[14] = pt_z + dt_z / 2;
							vertices[17] = mGrid.getDensity(glm::vec3(vertices[12], pt_y, vertices[14]));

							vertices[18] = pt_x + dt_x / 2;
							vertices[20] = pt_z + dt_z / 2;
							vertices[23] = mGrid.getDensity(glm::vec3(vertices[18], pt_y, vertices[20]));

							for (int k = 0; k <= 18; k += 6)
							{
								vertices[k] = (vertices[k] / mGrid.mD.mMax[2]);
								vertices[k + 1] = 0;
								vertices[k + 2] = (vertices[k + 2] / mGrid.mD.mMax[2]);
							}

							glBindVertexArray(VAO);
							glBindBuffer(GL_ARRAY_BUFFER, VBO);
							glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
							glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
							glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
							glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)0);
							glEnableVertexAttribArray(0);
							glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)(3 * sizeof(float)));
							glEnableVertexAttribArray(1);
							glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *)(5 * sizeof(float)));
							glEnableVertexAttribArray(2);
							glBindVertexArray(0);

							gridShader->use();
							glBindTexture(GL_TEXTURE_2D, smokeTexture);
							glUniform1i(glGetUniformLocation(gridShader->getId(), "mTexture"), 0);

							glm::mat4 view = Glb::Camera::getInstance().GetView();
							glm::mat4 projection = Glb::Camera::getInstance().GetProjection();

							glm::mat4 model = glm::mat4(1.0f);
							model = glm::translate(model, glm::vec3(0.0f, j, 0.0f));
							glUniformMatrix4fv(glGetUniformLocation(gridShader->getId(), "view"), 1, GL_FALSE, glm::value_ptr(view));
							glUniformMatrix4fv(glGetUniformLocation(gridShader->getId(), "projection"), 1, GL_FALSE, glm::value_ptr(projection));
							glUniformMatrix4fv(glGetUniformLocation(gridShader->getId(), "model"), 1, GL_FALSE, glm::value_ptr(model));
							gridShader->setFloat("contrast", Eulerian3dPara::contrast);
							glBindVertexArray(VAO);
							glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
						}
					}
				}
			}
		}

		GLuint Renderer::getTextureID()
		{
			return textureID;
		}

	}
}
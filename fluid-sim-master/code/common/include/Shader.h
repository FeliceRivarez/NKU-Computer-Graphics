#pragma once
#ifndef SHADER_H
#define SHADER_H

#include <glad/glad.h>
#include <string>
#include <glm/gtc/matrix_transform.hpp>

namespace Glb {
    class Shader {
    public:
        Shader();
        ~Shader();

        int32_t buildFromFile(std::string& vertPath, std::string& fragPath);
        int32_t buildFromFile(std::string& vertPath, std::string& fragPath, std::string& geomPath);
        void use();
        void unUse();
        GLuint getId();

        void setBool(const std::string& name, bool value);
        void setInt(const std::string& name, int value);
        void setFloat(const std::string& name, float value);
        void setVec2(const std::string& name, const glm::vec2& value);
        void setVec2(const std::string& name, float x, float y);
        void setVec3(const std::string& name, const glm::vec3& value);
        void setVec3(const std::string& name, float x, float y, float z);
        void setVec4(const std::string& name, const glm::vec4& value);
        void setVec4(const std::string& name, float x, float y, float z, float w);
        void setMat2(const std::string& name, const glm::mat2& mat);
        void setMat3(const std::string& name, const glm::mat3& mat);
        void setMat4(const std::string& name, const glm::mat4& mat);

    private:
        GLuint mId = 0; // shader program id

    };
}



#endif // !SHADER_H

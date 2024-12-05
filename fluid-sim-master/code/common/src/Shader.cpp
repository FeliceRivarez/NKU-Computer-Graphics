#include <Shader.h>
#include <fstream>
#include <sstream>
#include <iostream>


namespace Glb {
    Shader::Shader() {
        mId = 0;
    }

    Shader::~Shader() {
        if (mId != 0) {
            glDeleteProgram(mId);
        }
    }

    int32_t Shader::buildFromFile(std::string& vertPath, std::string& fragPath) {
        std::string vertexCode, fragmentCode;
        std::ifstream vertShaderFile, fragShaderFile;
        std::stringstream vShaderStream, fShaderStream;

        // read vertex shader
        vertShaderFile.open(vertPath);
        if (!vertShaderFile) {
            std::cout << "ERROR: Vertex shader file open failed" << std::endl;
            return -1;
        }
        vShaderStream << vertShaderFile.rdbuf();
        vertShaderFile.close();
        vertexCode = vShaderStream.str();
        const char* vShaderCode = vertexCode.c_str();   // saved

        // create vertex shader
        unsigned int vertex;
        vertex = glCreateShader(GL_VERTEX_SHADER);  // 着色器类型:GL_VERTEX_SHADER
        // set the shader's source code
        glShaderSource(vertex, 1, &vShaderCode, NULL);
        // compile the shader
        glCompileShader(vertex);

        // check compile errors
        int success;
        char infoLog[1024];
        glGetShaderiv(vertex, GL_COMPILE_STATUS, &success);
        if (!success)
        {
            glGetShaderInfoLog(vertex, 1024, NULL, infoLog);
            std::cout << "ERROR::VERTEX_SHADER_COMPILATION_ERROR: " << "\n" << infoLog << std::endl;
            return -1;
        }

        // read fragment Shader
        fragShaderFile.open(fragPath);
        if (!fragShaderFile) {
            std::cout << "ERROR: Fragment shader file open failed" << std::endl;
            return -1;
        }
        fShaderStream << fragShaderFile.rdbuf();
        fragShaderFile.close();
        fragmentCode = fShaderStream.str();
        const char* fShaderCode = fragmentCode.c_str();

        // create fragment Shader
        unsigned int fragment;
        fragment = glCreateShader(GL_FRAGMENT_SHADER);  // 着色器类型:GL_FRAGMENT_SHADER
        glShaderSource(fragment, 1, &fShaderCode, NULL);
        glCompileShader(fragment);

        // check compile errors
        glGetShaderiv(fragment, GL_COMPILE_STATUS, &success);
        if (!success)
        {
            glGetShaderInfoLog(fragment, 1024, NULL, infoLog);
            std::cout << "ERROR::FRAGMENT_SHADER_COMPILATION_ERROR: " << "\n" << infoLog << std::endl;
            return -1;
        }

        // create shader program
        mId = glCreateProgram();

        // attach shader to the program
        glAttachShader(mId, vertex);
        glAttachShader(mId, fragment);

        // link the shaders to be a executable program
        glLinkProgram(mId);

        // check link errors
        glGetProgramiv(mId, GL_LINK_STATUS, &success);
        if (!success)
        {
            glGetProgramInfoLog(mId, 1024, NULL, infoLog);
            std::cout << "ERROR::PROGRAM_LINKING_ERROR:" << "\n" << infoLog << std::endl;
            return -1;
        }

        // we already have the executable program, then just delete the shaders
        glDeleteShader(vertex);
        glDeleteShader(fragment);

        //std::cout << "success" << std::endl;

        return 0;
    }

    int32_t Shader::buildFromFile(std::string& vertPath, std::string& fragPath, std::string& geomPath) {
        std::string vertexCode, geometryCode, fragmentCode;
        std::ifstream vertShaderFile, geomShaderFile, fragShaderFile;
        std::stringstream vShaderStream, gShaderStream, fShaderStream;

        // read vertex shader
        vertShaderFile.open(vertPath);
        if (!vertShaderFile) {
            std::cout << "ERROR: Vertex shader file open failed" << std::endl;
            return -1;
        }
        vShaderStream << vertShaderFile.rdbuf();
        vertShaderFile.close();
        vertexCode = vShaderStream.str();
        const char* vShaderCode = vertexCode.c_str();   // saved

        // create vertex shader
        unsigned int vertex;
        vertex = glCreateShader(GL_VERTEX_SHADER);  // 着色器类型:GL_VERTEX_SHADER
        // set the shader's source code
        glShaderSource(vertex, 1, &vShaderCode, NULL);
        // compile the shader
        glCompileShader(vertex);

        // check compile errors
        int success;
        char infoLog[1024];
        glGetShaderiv(vertex, GL_COMPILE_STATUS, &success);
        if (!success)
        {
            glGetShaderInfoLog(vertex, 1024, NULL, infoLog);
            std::cout << "ERROR::VERTEX_SHADER_COMPILATION_ERROR: " << "\n" << infoLog << std::endl;
            return -1;
        }

        // read fragment Shader
        fragShaderFile.open(fragPath);
        if (!fragShaderFile) {
            std::cout << "ERROR: Fragment shader file open failed" << std::endl;
            return -1;
        }
        fShaderStream << fragShaderFile.rdbuf();
        fragShaderFile.close();
        fragmentCode = fShaderStream.str();
        const char* fShaderCode = fragmentCode.c_str();

        // create fragment Shader
        unsigned int fragment;
        fragment = glCreateShader(GL_FRAGMENT_SHADER);  // 着色器类型:GL_FRAGMENT_SHADER
        glShaderSource(fragment, 1, &fShaderCode, NULL);
        glCompileShader(fragment);
        
        // check compile errors
        glGetShaderiv(fragment, GL_COMPILE_STATUS, &success);
        if (!success)
        {
            glGetShaderInfoLog(fragment, 1024, NULL, infoLog);
            std::cout << "ERROR::FRAGMENT_SHADER_COMPILATION_ERROR: " << "\n" << infoLog << std::endl;
            return -1;
        }

        unsigned int geometry = 0;
        // if geometry shader exists
        if (geomPath != std::string("")) {

            // read geometry Shader
            geomShaderFile.open(geomPath);
            if (!geomShaderFile) {
                std::cout << "ERROR: Geometry shader file open failed" << std::endl;
                return -1;
            }
            gShaderStream << geomShaderFile.rdbuf();
            geomShaderFile.close();
            geometryCode = gShaderStream.str();
            const char* gShaderCode = geometryCode.c_str();
            
            // create geometry Shader
            geometry = glCreateShader(GL_GEOMETRY_SHADER);
            glShaderSource(geometry, 1, &gShaderCode, NULL);
            glCompileShader(geometry);
            
            // check compile errors
            glGetShaderiv(geometry, GL_COMPILE_STATUS, &success);
            if (!success)
            {
                glGetShaderInfoLog(geometry, 1024, NULL, infoLog);
                std::cout << "ERROR::GEOMETRY_SHADER_COMPILATION_ERROR: " << "\n" << infoLog << std::endl;
                return -1;
            }
        }

        // create shader program
        mId = glCreateProgram();

        // attach shader to the program
        glAttachShader(mId, vertex);
        glAttachShader(mId, fragment);
        if (geomPath != std::string("")) {
            glAttachShader(mId, geometry);
        }

        // link the shaders to be a executable program
        glLinkProgram(mId);

        // check link errors
        glGetProgramiv(mId, GL_LINK_STATUS, &success);
        if (!success)
        {
            glGetProgramInfoLog(mId, 1024, NULL, infoLog);
            std::cout << "ERROR::PROGRAM_LINKING_ERROR:" << "\n" << infoLog << std::endl;
            return -1;
        }

        // we already have the executable program, then just delete the shaders
        glDeleteShader(vertex);
        glDeleteShader(fragment);
        if (geomPath != std::string("")) {
            glDeleteShader(geometry);
        }
        std::cout << "success" << std::endl;

        return 0;
    }

    void Shader::use() {
        glUseProgram(mId);
    }

    void Shader::unUse() {
        glUseProgram(0);
    }

    GLuint Shader::getId() {
        return mId;
    }

    void Shader::setBool(const std::string& name, bool value)
    {
        glUniform1i(glGetUniformLocation(mId, name.c_str()), (int)value);
    }

    void Shader::setInt(const std::string& name, int value)
    {
        glUniform1i(glGetUniformLocation(mId, name.c_str()), value);
    }

    void Shader::setFloat(const std::string& name, float value)
    {
        glUniform1f(glGetUniformLocation(mId, name.c_str()), value);
    }

    void Shader::setVec2(const std::string& name, const glm::vec2& value)
    {
        glUniform2fv(glGetUniformLocation(mId, name.c_str()), 1, &value[0]);
    }
    void Shader::setVec2(const std::string& name, float x, float y)
    {
        glUniform2f(glGetUniformLocation(mId, name.c_str()), x, y);
    }

    void Shader::setVec3(const std::string& name, const glm::vec3& value)
    {
        glUniform3fv(glGetUniformLocation(mId, name.c_str()), 1, &value[0]);
    }
    void Shader::setVec3(const std::string& name, float x, float y, float z)
    {
        glUniform3f(glGetUniformLocation(mId, name.c_str()), x, y, z);
    }

    void Shader::setVec4(const std::string& name, const glm::vec4& value)
    {
        glUniform4fv(glGetUniformLocation(mId, name.c_str()), 1, &value[0]);
    }
    void Shader::setVec4(const std::string& name, float x, float y, float z, float w)
    {
        glUniform4f(glGetUniformLocation(mId, name.c_str()), x, y, z, w);
    }

    void Shader::setMat2(const std::string& name, const glm::mat2& mat)
    {
        glUniformMatrix2fv(glGetUniformLocation(mId, name.c_str()), 1, GL_FALSE, &mat[0][0]);
    }

    void Shader::setMat3(const std::string& name, const glm::mat3& mat)
    {
        glUniformMatrix3fv(glGetUniformLocation(mId, name.c_str()), 1, GL_FALSE, &mat[0][0]);
    }

    void Shader::setMat4(const std::string& name, const glm::mat4& mat)
    {
        glUniformMatrix4fv(glGetUniformLocation(mId, name.c_str()), 1, GL_FALSE, &mat[0][0]);
    }
}



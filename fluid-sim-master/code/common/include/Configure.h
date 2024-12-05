#pragma once
#ifndef __CONFIGURE_H__
#define __CONFIGURE_H__

#include <iostream>
#include <string>
#include "Component.h"
#include <vector>
#include "glm/glm.hpp"

#define LERP(a, b, t) (1 - t) * a + t *b

#ifndef __MINMAX_DEFINED
#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))
#endif

extern int imageWidth;
extern int imageHeight;

extern int windowWidth;
extern int windowHeight;

extern float fontSize;

extern bool simulating;

namespace Eulerian2dPara
{
    struct SourceSmoke {
        glm::ivec2 position = glm::ivec2(0);
        glm::vec2 velocity = glm::vec2(0.0f);
        float density = 0.0f;
        float temp = 0.0f;
    };

    extern int theDim2d[];
    extern std::vector<SourceSmoke> source;
    extern float theCellSize2d;
    extern bool addSolid;

    extern float dt;

    extern float contrast;
    extern int drawModel;
    extern int gridNum;


    extern float airDensity;
    extern float ambientTemp;
    extern float boussinesqAlpha;
    extern float boussinesqBeta;
    extern float vorticityConst;
}

namespace Eulerian3dPara
{
    struct SourceSmoke {
        glm::ivec3 position = glm::ivec3(0);
        glm::vec3 velocity = glm::vec3(0.0f);
        float density = 0.0f;
        float temp = 0.0f;
    };

    extern int theDim3d[];
    extern float theCellSize3d;
    extern std::vector<SourceSmoke> source;
    extern bool addSolid;

    extern float contrast;
    extern bool oneSheet;
    extern float distanceX;
    extern float distanceY;
    extern float distanceZ;
    extern bool xySheetsON;
    extern bool yzSheetsON;
    extern bool xzSheetsON;
    extern int xySheetsNum;
    extern int yzSheetsNum;
    extern int xzSheetsNum;
    extern int drawModel;
    extern int gridNumX;
    extern int gridNumY;
    extern int gridNumZ;

    extern float dt;

    extern float airDensity;
    extern float ambientTemp;
    extern float boussinesqAlpha;
    extern float boussinesqBeta;
    extern float vorticityConst;

}

namespace Lagrangian2dPara
{
    struct FluidBlock {
        glm::vec2 lowerCorner = glm::vec2(0.0f, 0.0f);
        glm::vec2 upperCorner = glm::vec2(0.0f, 0.0f);
        glm::vec2 initVel = glm::vec2(0.0f, 0.0f);
        float particleSpace = 0.02f;
    };

    extern float scale;
    extern std::vector<FluidBlock> fluidBlocks;

    extern float dt;
    extern int substep;
    extern float maxVelocity;
    extern float velocityAttenuation;
    extern float eps;

    extern float supportRadius;
    extern float particleRadius;
    extern float particleDiameter;
    extern float gravityX;
    extern float gravityY;
    extern float density;
    extern float stiffness;
    extern float exponent;
    extern float viscosity;
}

namespace Lagrangian3dPara
{

    struct FluidBlock {
        glm::vec3 lowerCorner = glm::vec3(0.0f, 0.0f, 0.0f);
        glm::vec3 upperCorner = glm::vec3(0.0f, 0.0f, 0.0f);
        glm::vec3 initVel = glm::vec3(0.0f, 0.0f, 0.0f);
        float particleSpace = 0.02f;
    };

    extern float scale;
    extern std::vector<FluidBlock> fluidBlocks;


    extern float dt;
    extern int substep;
    extern float maxVelocity;
    extern float velocityAttenuation;
    extern float eps;

    extern float supportRadius;
    extern float particleRadius;
    extern float particleDiameter;

    extern float gravityX;
    extern float gravityY;
    extern float gravityZ;

    extern float density;
    extern float stiffness;
    extern float exponent;
    extern float viscosity;
}

extern std::string shaderPath;
extern std::string picturePath;

extern std::vector<Glb::Component *> methodComponents;

#endif // !__CONFIGURE_H__

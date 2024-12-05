#pragma once
#ifndef __MACGRID_2D_H__
#define __MACGRID_2D_H__

#include <windows.h>
#include <glm/glm.hpp>
#include "GridData2d.h"
#include <Logger.h>

namespace FluidSimulation
{
    namespace Eulerian2d
    {

        class MACGrid2d
        {
        public:
            MACGrid2d();
            ~MACGrid2d();
            MACGrid2d(const MACGrid2d &orig);
            MACGrid2d &operator=(const MACGrid2d &orig);

            void reset();

            glm::vec4 getRenderColor(int i, int j);
            glm::vec4 getRenderColor(const glm::vec2 &pt);

            // Setup
            void initialize();
            void createSolids();
            void updateSources();

            // advect
            glm::vec2 semiLagrangian(const glm::vec2 &pt, double dt);

            // get value
            glm::vec2 getVelocity(const glm::vec2 &pt);
            double getVelocityX(const glm::vec2 &pt);
            double getVelocityY(const glm::vec2 &pt);
            double getTemperature(const glm::vec2 &pt);
            double getDensity(const glm::vec2 &pt);

            enum Direction
            {
                X,
                Y
            };

            // get point
            glm::vec2 getCenter(int i, int j);
            glm::vec2 getLeft(int i, int j);
            glm::vec2 getRight(int i, int j);
            glm::vec2 getTop(int i, int j);
            glm::vec2 getBottom(int i, int j);

            void getCell(int index, int &i, int &j);
            int getIndex(int i, int j);
            bool isNeighbor(int i0, int j0, int i1, int j1);
            bool isValid(int i, int j, Direction d);

            int isSolidCell(int i, int j);              // Returns 1 if true, else otherwise
            int isSolidFace(int i, int j, Direction d); // Returns 1 if true, else otherwise

            bool inSolid(const glm::vec2 &pt);
            bool inSolid(const glm::vec2 &pt, int &i, int &j);

            bool intersects(const glm::vec2 &pt, const glm::vec2 &dir, int i, int j, double &time);
            int numSolidCells();

            // 计算两个网格单元之间的压力系数
            double getPressureCoeffBetweenCells(int i0, int j0, int i1, int j1);
            
            // 计算散度
            double getDivergence(int i, int j);
            // 检查散度
            double checkDivergence(int i, int j);
            bool checkDivergence();

            // 计算Boussinesq Force
            double getBoussinesqForce(const glm::vec2 &pt);

            float cellSize;
            int dim[2];

            Glb::GridData2dX mU;        // x velocity
            Glb::GridData2dY mV;        // y velocity
            Glb::CubicGridData2d mD;    // density
            Glb::CubicGridData2d mT;    // temperature

            Glb::GridData2d mSolid;     // solid
        };

#define FOR_EACH_CELL                                                \
    for (int j = 0; j < Eulerian2dPara::theDim2d[MACGrid2d::Y]; j++) \
        for (int i = 0; i < Eulerian2dPara::theDim2d[MACGrid2d::X]; i++)

#define FOR_EACH_LINE                                                    \
    for (int j = 0; j < Eulerian2dPara::theDim2d[MACGrid2d::Y] + 1; j++) \
        for (int i = 0; i < Eulerian2dPara::theDim2d[MACGrid2d::X] + 1; i++)

    }
}

#endif
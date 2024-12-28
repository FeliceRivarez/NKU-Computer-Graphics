#include "MACGrid2d.h"
#include "Configure.h"
#include <math.h>
#include <map>
#include <stdio.h>

namespace FluidSimulation
{
    namespace Eulerian2d
    {

        MACGrid2d::MACGrid2d()
        {
            cellSize = Eulerian2dPara::theCellSize2d;
            dim[0] = Eulerian2dPara::theDim2d[0];
            dim[1] = Eulerian2dPara::theDim2d[1];
            initialize();
        }

        MACGrid2d::MACGrid2d(const MACGrid2d &orig)
        {
            mU = orig.mU;
            mV = orig.mV;
            mD = orig.mD;
            mT = orig.mT;
            mSolid = orig.mSolid;
        }

        MACGrid2d &MACGrid2d::operator=(const MACGrid2d &orig)
        {
            if (&orig == this)
            {
                return *this;
            }
            mU = orig.mU;
            mV = orig.mV;
            mD = orig.mD;
            mT = orig.mT;
            mSolid = orig.mSolid;

            return *this;
        }

        MACGrid2d::~MACGrid2d()
        {
        }

        void MACGrid2d::reset()
        {
            mU.initialize(0.0);
            mV.initialize(0.0);
            mD.initialize(0.0);
            mT.initialize(Eulerian2dPara::ambientTemp);
        }

        // void MACGrid2d::createSolids()
        // {
        //     mSolid.initialize();
        //     if (Eulerian2dPara::addSolid) {
        //         int j = dim[1] / 2;
        //         for (int i = dim[0] / 4; i < dim[0] * 3 / 4; i++) {
        //             mSolid(i, j) = 1;
        //         }
        //     }
        // }
        void MACGrid2d::createSolids()
        {
            mSolid.initialize();
            if (Eulerian2dPara::addSolid) {
                int centerX = dim[0] / 2;  // 圆心 x 坐标
                int centerY = dim[1] / 2;  // 圆心 y 坐标
                int R = dim[0] / 4;   // 圆的半径

                for (int i = 0; i < dim[0]; i++) {
                    for (int j = 0; j < dim[1]; j++) {
                        // 计算点到圆心的距离平方
                        int len = (i - centerX) * (i - centerX) + (j - centerY) * (j - centerY);
                        int out = R * R;
                        int inn = (R - 1) * (R - 1);

                        // 判断点是否在圆的边界上
                        if (len <= out && len >= inn) {
                            // 排除底部和顶部开口区域
                            if (j > centerY - R + 2 && j<centerY+R-7 ||j>centerY+R-5) {
                                mSolid(i, j) = 1;
                            }
                        }
                    }
                }
            }
        }

        void MACGrid2d::updateSources()
        {
            for (int i = 0; i < Eulerian2dPara::source.size(); i++) {
                int x = Eulerian2dPara::source[i].position.x;
                int y = Eulerian2dPara::source[i].position.y;
                mT(x, y) = Eulerian2dPara::source[i].temp;
                mD(x, y) = Eulerian2dPara::source[i].density;
                mU(x, y) = Eulerian2dPara::source[i].velocity.x;
                mV(x, y) = Eulerian2dPara::source[i].velocity.y;
            }
        }

        void MACGrid2d::initialize()
        {
            reset();
            createSolids();
            assert(checkDivergence());
        }

        // Boussinesq Force
        double MACGrid2d::getBoussinesqForce(const glm::vec2 &pos)
        {
            double temperature = getTemperature(pos);
            double smokeDensity = getDensity(pos);

            double yforce = -Eulerian2dPara::boussinesqAlpha * smokeDensity +
                            Eulerian2dPara::boussinesqBeta * (temperature - Eulerian2dPara::ambientTemp);

            return yforce;
        }


        // 计算散度
        double MACGrid2d::getDivergence(int i, int j)
        {

            double x1 = isSolidCell(i + 1, j) ? 0.0 : mU(i + 1, j);
            double x0 = isSolidCell(i - 1, j) ? 0.0 : mU(i, j);

            double y1 = isSolidCell(i, j + 1) ? 0.0 : mV(i, j + 1);
            double y0 = isSolidCell(i, j - 1) ? 0.0 : mV(i, j);

            double xdiv = x1 - x0;
            double ydiv = y1 - y0;
            double div = (xdiv + ydiv) / cellSize;

            return div;
        }

        double MACGrid2d::checkDivergence(int i, int j)
        {
            double x1 = mU(i + 1, j);
            double x0 = mU(i, j);
            double y1 = mV(i, j + 1);
            double y0 = mV(i, j);
            double xdiv = x1 - x0;
            double ydiv = y1 - y0;
            double div = (xdiv + ydiv) / cellSize;
            return div;
        }

        bool MACGrid2d::checkDivergence()
        {
            FOR_EACH_CELL
            {
                double div = checkDivergence(i, j);
                if (fabs(div) > 0.01)
                {
                    return false;
                }
            }
            return true;
        }

        glm::vec2 MACGrid2d::semiLagrangian(const glm::vec2 &pt, double dt)
        {
            glm::vec2 vel = getVelocity(pt);
            glm::vec2 pos = pt - vel * (float)dt;

            pos[0] = max(0.0, min((dim[0] - 1) * cellSize, pos[0]));
            pos[1] = max(0.0, min((dim[1] - 1) * cellSize, pos[1]));

            int i, j;
            if (inSolid(pt, i, j))
            {
                double t = 0;
                if (intersects(pt, vel, i, j, t))
                {
                    pos = pt - vel * (float)t;
                }
                else
                {
                    Glb::Logger::getInstance().addLog("Error: something goes wrong during advection");
                }
            }
            return pos;
        }

        // 获取相交的具体时间
        bool MACGrid2d::intersects(const glm::vec2 &wPos, const glm::vec2 &wDir, int i, int j, double &time)
        {
            glm::vec2 pos = getCenter(i, j);

            glm::vec2 rayStart = wPos - pos;
            glm::vec2 rayDir = wDir;

            double tmin = -9999999999.0;
            double tmax = 9999999999.0;

            double min = -0.5 * cellSize;
            double max = 0.5 * cellSize;

            for (int i = 0; i < 3; i++)
            {
                double e = rayStart[i];
                double f = rayDir[i];
                if (fabs(f) > 0.000000001)
                {
                    double t1 = (min - e) / f;
                    double t2 = (max - e) / f;
                    if (t1 > t2)
                        std::swap(t1, t2);
                    if (t1 > tmin)
                        tmin = t1;
                    if (t2 < tmax)
                        tmax = t2;
                    if (tmin > tmax)
                        return false;
                    if (tmax < 0)
                        return false;
                }
                else if (e < min || e > max)
                    return false;
            }
            if (tmin >= 0)
            {
                time = tmin;
                return true;
            }
            else
            {
                time = tmax;
                return true;
            }
            return false;
        }


        int MACGrid2d::getIndex(int i, int j)
        {
            if (i < 0 || i > dim[0] - 1)
                return -1;
            if (j < 0 || j > dim[1] - 1)
                return -1;

            int col = i;
            int row = j * dim[0];
            return col + row;
        }


        void MACGrid2d::getCell(int index, int &i, int &j)
        {
            j = (int)index / dim[0]; // row
            i = index - j * dim[0];  // col
        }


        glm::vec2 MACGrid2d::getCenter(int i, int j)
        {
            double xstart = cellSize / 2.0;
            double ystart = cellSize / 2.0;

            double x = xstart + i * cellSize;
            double y = ystart + j * cellSize;
            return glm::vec2(x, y);
        }

        glm::vec2 MACGrid2d::getLeft(int i, int j)
        {
            return getCenter(i, j) - glm::vec2(cellSize * 0.5, 0.0);
        }

        glm::vec2 MACGrid2d::getRight(int i, int j)
        {
            return getCenter(i, j) + glm::vec2(cellSize * 0.5, 0.0);
        }

        glm::vec2 MACGrid2d::getTop(int i, int j)
        {
            return getCenter(i, j) + glm::vec2(0.0, cellSize * 0.5);
        }

        glm::vec2 MACGrid2d::getBottom(int i, int j)
        {
            return getCenter(i, j) - glm::vec2(0.0, cellSize * 0.5);
        }

        glm::vec2 MACGrid2d::getVelocity(const glm::vec2 &pt)
        {
            if (inSolid(pt))
            {
                return glm::vec2(0, 0);
            }

            glm::vec2 vel;
            vel[0] = getVelocityX(pt);
            vel[1] = getVelocityY(pt);
            return vel;
        }

        double MACGrid2d::getVelocityX(const glm::vec2 &pt)
        {
            return mU.interpolate(pt);
        }

        double MACGrid2d::getVelocityY(const glm::vec2 &pt)
        {
            return mV.interpolate(pt);
        }

        double MACGrid2d::getTemperature(const glm::vec2 &pt)
        {
            return mT.interpolate(pt);
        }

        double MACGrid2d::getDensity(const glm::vec2 &pt)
        {
            return mD.interpolate(pt);
        }

        int MACGrid2d::numSolidCells()
        {
            int numSolid = 0;
            FOR_EACH_CELL { numSolid += mSolid(i, j); }
            return numSolid;
        }

        bool MACGrid2d::inSolid(const glm::vec2 &pt)
        {
            int i, j, k;
            mSolid.getCell(pt, i, j);
            return isSolidCell(i, j) == 1;
        }

        bool MACGrid2d::inSolid(const glm::vec2 &pt, int &i, int &j)
        {
            mSolid.getCell(pt, i, j);
            return isSolidCell(i, j) == 1;
        }

        int MACGrid2d::isSolidCell(int i, int j)
        {
            bool containerBoundary = (i < 0 || i > dim[0] - 1) ||
                                     (j < 0 || j > dim[1] - 1);

            bool objectBoundary = (mSolid(i, j) == 1);

            return containerBoundary || objectBoundary ? 1 : 0;
        }

        int MACGrid2d::isSolidFace(int i, int j, MACGrid2d::Direction d)
        {
            if (d == X && (i == 0 || i == dim[0]))
                return 1;
            else if (d == Y && (j == 0 || j == dim[1]))
                return 1;

            if (d == X && (mSolid(i, j) || mSolid(i - 1, j)))
                return 1;
            if (d == Y && (mSolid(i, j) || mSolid(i, j - 1)))
                return 1;

            return 0;
        }


        bool MACGrid2d::isNeighbor(int i0, int j0, int i1, int j1)
        {
            if (abs(i0 - i1) == 1 && j0 == j1)
                return true;
            if (abs(j0 - j1) == 1 && i0 == i1)
                return true;
            return false;
        }


        double MACGrid2d::getPressureCoeffBetweenCells(
            int i, int j, int pi, int pj)
        {
            // 同一个cell
            if (i == pi && j == pj) // self
            {
                int numSolidNeighbors = (isSolidCell(i + 1, j) +
                                         isSolidCell(i - 1, j) +
                                         isSolidCell(i, j + 1) +
                                         isSolidCell(i, j - 1));
                // Return number of non-solid boundaries around cel ij
                return 4.0 - numSolidNeighbors; // 二维，4
            }
            // 如果两个是邻居且后者为固体，error
            if (isNeighbor(i, j, pi, pj) && !isSolidCell(pi, pj))
                return -1.0;
            // 不是邻居
            return 0.0;
        }


        glm::vec4 MACGrid2d::getRenderColor(int i, int j)
        {
            double value = mD(i, j);
            return glm::vec4(1.0, 1.0, 1.0, value);
        }


        glm::vec4 MACGrid2d::getRenderColor(const glm::vec2 &pt)
        {
            double value = getDensity(pt);
            return glm::vec4(value, value, value, value);
        }


        // 确保在界内
        bool MACGrid2d::isValid(int i, int j, MACGrid2d::Direction d)
        {
            switch (d)
            {
            case X:
                return (i >= 0 && i < dim[X] + 1 &&
                        j >= 0 && j < dim[Y]);
            case Y:
                return (i >= 0 && i < dim[X] &&
                        j >= 0 && j < dim[Y] + 1);
            }
            Glb::Logger::getInstance().addLog("Error: bad direction");
            return false;
        }
    }
}
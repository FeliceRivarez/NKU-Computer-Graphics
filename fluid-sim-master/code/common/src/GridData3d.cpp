#include "GridData3d.h"
#include "Configure.h"

namespace Glb
{

    GridData3d::GridData3d() : mDfltValue(0.0), mMax(0.0, 0.0, 0.0), cellSize(Eulerian3dPara::theCellSize3d)
    {
        dim[0] = Eulerian3dPara::theDim3d[0];
        dim[1] = Eulerian3dPara::theDim3d[1];
        dim[2] = Eulerian3dPara::theDim3d[2];
    }

    GridData3d::GridData3d(const GridData3d &orig) : mDfltValue(orig.mDfltValue)
    {
        mData = orig.mData;
        mMax = orig.mMax;
        cellSize = orig.cellSize;
        dim[0] = orig.dim[0];
        dim[1] = orig.dim[1];
        dim[2] = orig.dim[2];
    }

    GridData3d::~GridData3d()
    {
    }

    ublas::vector<double> &GridData3d::data()
    {
        return mData;
    }

    GridData3d &GridData3d::operator=(const GridData3d &orig)
    {
        if (this == &orig)
        {
            return *this;
        }
        mData = orig.mData;
        mMax = orig.mMax;
        cellSize = orig.cellSize;
        dim[0] = orig.dim[0];
        dim[1] = orig.dim[1];
        dim[2] = orig.dim[2];
        return *this;
    }

    void GridData3d::initialize(double dfltValue)
    {
        mDfltValue = dfltValue;
        mMax[0] = cellSize * dim[0];
        mMax[1] = cellSize * dim[1];
        mMax[2] = cellSize * dim[2];
        mData.resize(dim[0] * dim[1] * dim[2], false);
        std::fill(mData.begin(), mData.end(), mDfltValue);
    }

    double &GridData3d::operator()(int i, int j, int k)
    {
        static double dflt = 0;
        dflt = mDfltValue;

        if (i < 0 || j < 0 || k < 0 ||
            i > dim[0] - 1 ||
            j > dim[1] - 1 ||
            k > dim[2] - 1)
            return dflt;

        int col = i;
        int row = k * dim[0];
        int stack = j * dim[0] * dim[2];

        return mData(col + row + stack);
    }

    void GridData3d::getCell(const glm::vec3 &pt, int &i, int &j, int &k)
    {
        glm::vec3 pos = worldToSelf(pt);
        i = (int)(pos[0] / cellSize);
        j = (int)(pos[1] / cellSize);
        k = (int)(pos[2] / cellSize);
    }

    double GridData3d::interpolate(const glm::vec3 &pt)
    {
        glm::vec3 pos = worldToSelf(pt);

        int i = (int)(pos[0] / cellSize);
        int j = (int)(pos[1] / cellSize);
        int k = (int)(pos[2] / cellSize);

        double scale = 1.0 / cellSize;
        double fractx = scale * (pos[0] - i * cellSize);
        double fracty = scale * (pos[1] - j * cellSize);
        double fractz = scale * (pos[2] - k * cellSize);

        assert(fractx < 1.0 && fractx >= 0);
        assert(fracty < 1.0 && fracty >= 0);
        assert(fractz < 1.0 && fractz >= 0);

        double tmp1 = (*this)(i, j, k);
        double tmp2 = (*this)(i, j + 1, k);
        double tmp3 = (*this)(i + 1, j, k);
        double tmp4 = (*this)(i + 1, j + 1, k);

        double tmp5 = (*this)(i, j, k + 1);
        double tmp6 = (*this)(i, j + 1, k + 1);
        double tmp7 = (*this)(i + 1, j, k + 1);
        double tmp8 = (*this)(i + 1, j + 1, k + 1);

        double tmp12 = LERP(tmp1, tmp2, fracty);
        double tmp34 = LERP(tmp3, tmp4, fracty);

        double tmp56 = LERP(tmp5, tmp6, fracty);
        double tmp78 = LERP(tmp7, tmp8, fracty);

        double tmp1234 = LERP(tmp12, tmp34, fractx);
        double tmp5678 = LERP(tmp56, tmp78, fractx);

        double tmp = LERP(tmp1234, tmp5678, fractz);
        return tmp;
    }

    glm::vec3 GridData3d::worldToSelf(const glm::vec3 &pt) const
    {
        glm::vec3 out;
        out[0] = min(max(0.0, pt[0] - cellSize * 0.5), mMax[0]);
        out[1] = min(max(0.0, pt[1] - cellSize * 0.5), mMax[1]);
        out[2] = min(max(0.0, pt[2] - cellSize * 0.5), mMax[2]);
        return out;
    }

    GridData3dX::GridData3dX() : GridData3d()
    {
    }

    GridData3dX::~GridData3dX()
    {
    }

    void GridData3dX::initialize(double dfltValue)
    {
        GridData3d::initialize(dfltValue);
        mMax[0] = cellSize * (dim[0] + 1);
        mMax[1] = cellSize * dim[1];
        mMax[2] = cellSize * dim[2];
        mData.resize((dim[0] + 1) * dim[1] * dim[2], false);
        std::fill(mData.begin(), mData.end(), mDfltValue);
    }

    double &GridData3dX::operator()(int i, int j, int k)
    {
        static double dflt = 0;
        dflt = mDfltValue; // Protect against setting the default value

        if (i < 0 || i > dim[0])
            return dflt;

        if (j < 0)
            j = 0;
        if (j > dim[1] - 1)
            j = dim[1] - 1;
        if (k < 0)
            k = 0;
        if (k > dim[2] - 1)
            k = dim[2] - 1;

        int col = i;
        int row = k * (dim[0] + 1);
        int stack = j * (dim[0] + 1) * dim[2];
        return mData(stack + row + col);
    }

    glm::vec3 GridData3dX::worldToSelf(const glm::vec3 &pt) const
    {
        glm::vec3 out;
        out[0] = min(max(0.0, pt[0]), mMax[0]);
        out[1] = min(max(0.0, pt[1] - cellSize * 0.5), mMax[1]);
        out[2] = min(max(0.0, pt[2] - cellSize * 0.5), mMax[2]);
        return out;
    }

    GridData3dY::GridData3dY() : GridData3d()
    {
    }

    GridData3dY::~GridData3dY()
    {
    }

    void GridData3dY::initialize(double dfltValue)
    {
        GridData3d::initialize(dfltValue);
        mMax[0] = cellSize * dim[0];
        mMax[1] = cellSize * (dim[1] + 1);
        mMax[2] = cellSize * dim[2];
        mData.resize(dim[0] * (dim[1] + 1) * dim[2], false);
        std::fill(mData.begin(), mData.end(), mDfltValue);
    }

    double &GridData3dY::operator()(int i, int j, int k)
    {
        static double dflt = 0;
        dflt = mDfltValue;

        if (j < 0 || j > dim[1])
            return dflt;

        if (i < 0)
            i = 0;
        if (i > dim[0] - 1)
            i = dim[0] - 1;
        if (k < 0)
            k = 0;
        if (k > dim[2] - 1)
            k = dim[2] - 1;

        int col = i;
        int row = k * dim[0];
        int stack = j * dim[0] * dim[2];
        return mData(stack + row + col);
    }

    glm::vec3 GridData3dY::worldToSelf(const glm::vec3 &pt) const
    {
        glm::vec3 out;
        out[0] = min(max(0.0, pt[0] - cellSize * 0.5), mMax[0]);
        out[1] = min(max(0.0, pt[1]), mMax[1]);
        out[2] = min(max(0.0, pt[2] - cellSize * 0.5), mMax[2]);
        return out;
    }

    GridData3dZ::GridData3dZ() : GridData3d()
    {
    }

    GridData3dZ::~GridData3dZ()
    {
    }

    void GridData3dZ::initialize(double dfltValue)
    {
        GridData3d::initialize(dfltValue);
        mMax[0] = cellSize * dim[0];
        mMax[1] = cellSize * dim[1];
        mMax[2] = cellSize * (dim[2] + 1);
        mData.resize(dim[0] * dim[1] * (dim[2] + 1), false);
        std::fill(mData.begin(), mData.end(), mDfltValue);
    }

    double &GridData3dZ::operator()(int i, int j, int k)
    {
        static double dflt = 0;
        dflt = mDfltValue;

        if (k < 0 || k > dim[2])
            return dflt;

        if (i < 0)
            i = 0;
        if (i > dim[0] - 1)
            i = dim[0] - 1;
        if (j < 0)
            j = 0;
        if (j > dim[1] - 1)
            j = dim[1] - 1;

        int col = i;
        int row = k * dim[0];
        int stack = j * dim[0] * (dim[2] + 1);

        return mData(stack + row + col);
    }

    glm::vec3 GridData3dZ::worldToSelf(const glm::vec3 &pt) const
    {
        glm::vec3 out;
        out[0] = min(max(0.0, pt[0] - cellSize * 0.5), mMax[0]);
        out[1] = min(max(0.0, pt[1] - cellSize * 0.5), mMax[1]);
        out[2] = min(max(0.0, pt[2]), mMax[2]);
        return out;
    }

    CubicGridData3d::CubicGridData3d() : GridData3d()
    {
    }

    CubicGridData3d::CubicGridData3d(const CubicGridData3d &orig) : GridData3d(orig)
    {
    }

    CubicGridData3d::~CubicGridData3d()
    {
    }

    double CubicGridData3d::cubic(double q1, double q2, double q3, double q4, double t)
    {
        double deltaq = q3 - q2;
        double d1 = (q3 - q1) * 0.5;
        double d2 = (q4 - q2) * 0.5;

        if (deltaq > 0.0001)
        {
            d1 = d1 > 0 ? d1 : 0.0;
            d2 = d2 > 0 ? d2 : 0.0;
        }
        else if (deltaq < 0.0001)
        {
            d1 = d1 < 0 ? d1 : 0.0;
            d2 = d2 < 0 ? d2 : 0.0;
        }

        double tmp = q2 + d1 * t + (3 * deltaq - 2 * d1 - d2) * t * t + (-2 * deltaq + d1 + d2) * t * t * t;
        return tmp;
    }

    double CubicGridData3d::interpY(int i, int j, int k, double fracty)
    {
        double tmp1 = (*this)(i, j - 1 < 0 ? j : j - 1, k);
        double tmp2 = (*this)(i, j, k);
        double tmp3 = (*this)(i, j + 1, k);
        double tmp4 = (*this)(i, j + 2, k);
        return cubic(tmp1, tmp2, tmp3, tmp4, fracty);
    }

    double CubicGridData3d::interpX(int i, int j, int k, double fracty, double fractx)
    {
        double tmp1 = interpY(i - 1 < 0 ? i : i - 1, j, k, fracty);
        double tmp2 = interpY(i, j, k, fracty);
        double tmp3 = interpY(i + 1, j, k, fracty);
        double tmp4 = interpY(i + 2, j, k, fracty);
        return cubic(tmp1, tmp2, tmp3, tmp4, fractx);
    }

    double CubicGridData3d::interpolate(const glm::vec3 &pt)
    {
        glm::vec3 pos = worldToSelf(pt);

        int i = (int)(pos[0] / cellSize);
        int j = (int)(pos[1] / cellSize);
        int k = (int)(pos[2] / cellSize);

        double scale = 1.0 / cellSize;
        double fractx = scale * (pos[0] - i * cellSize);
        double fracty = scale * (pos[1] - j * cellSize);
        double fractz = scale * (pos[2] - k * cellSize);

        assert(fractx < 1.0 && fractx >= 0);
        assert(fracty < 1.0 && fracty >= 0);
        assert(fractz < 1.0 && fractz >= 0);

        double tmp1 = interpX(i, j, k - 1 < 0 ? k : k - 1, fracty, fractx);
        double tmp2 = interpX(i, j, k, fracty, fractx);
        double tmp3 = interpX(i, j, k + 1, fracty, fractx);
        double tmp4 = interpX(i, j, k + 2, fracty, fractx);

        double tmp = cubic(tmp1, tmp2, tmp3, tmp4, fractz);
        return tmp;
    }
}

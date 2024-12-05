#pragma once
#ifndef __GRID_DATA_3D_H__
#define __GRID_DATA_3D_H__

#pragma warning(disable: 4244 4267 4996)
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <glm/glm.hpp>

namespace Glb {

	using namespace boost::numeric;

	class GridData3d
	{
	public:
		GridData3d();
		GridData3d(const GridData3d& orig);
		virtual ~GridData3d();
		virtual GridData3d& operator=(const GridData3d& orig);

		virtual void initialize(double dfltValue = 0.0);

		virtual double& operator()(int i, int j, int k);

		virtual double interpolate(const glm::vec3& pt);

		ublas::vector<double>& data();


		virtual void getCell(const glm::vec3& pt, int& i, int& j, int& k);

		glm::vec3 mMax;

	protected:

		virtual glm::vec3 worldToSelf(const glm::vec3& pt) const;
		double mDfltValue;
		ublas::vector<double> mData;
		float cellSize;
		int dim[3];
	};

	class GridData3dX : public GridData3d
	{
	public:
		GridData3dX();
		virtual ~GridData3dX();
		virtual void initialize(double dfltValue = 0.0);
		virtual double& operator()(int i, int j, int k);
		virtual glm::vec3 worldToSelf(const glm::vec3& pt) const;
	};

	class GridData3dY : public GridData3d
	{
	public:
		GridData3dY();
		virtual ~GridData3dY();
		virtual void initialize(double dfltValue = 0.0);
		virtual double& operator()(int i, int j, int k);
		virtual glm::vec3 worldToSelf(const glm::vec3& pt) const;
	};

	class GridData3dZ : public GridData3d
	{
	public:
		GridData3dZ();
		virtual ~GridData3dZ();
		virtual void initialize(double dfltValue = 0.0);
		virtual double& operator()(int i, int j, int k);
		virtual glm::vec3 worldToSelf(const glm::vec3& pt) const;
	};

	class CubicGridData3d : public GridData3d
	{
	public:
		CubicGridData3d();
		CubicGridData3d(const CubicGridData3d& orig);
		virtual ~CubicGridData3d();
		virtual double interpolate(const glm::vec3& pt);

	protected:

		double cubic(double q1, double q2, double q3, double q4, double t);
		double interpX(int i, int j, int k, double fracty, double fractx);
		double interpY(int i, int j, int k, double fracty);
	};
}

#endif
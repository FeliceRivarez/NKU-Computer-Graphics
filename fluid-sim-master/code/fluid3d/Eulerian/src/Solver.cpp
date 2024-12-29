#include "fluid3d/Eulerian/include/Solver.h"
#include "Configure.h"
#include "Global.h"

namespace FluidSimulation
{
    namespace Eulerian3d
    {
        Solver::Solver(MACGrid3d &grid) : mGrid(grid)
        {
            mGrid.reset();
        }

         void Solver::advection(MACGrid3d& grid){
            auto compute_X = grid.mU;
            auto compute_Y = grid.mV;
            auto compute_Z = grid.mW;
            auto compute_D = grid.mD;
            auto compute_T = grid.mT;

            //遍历网格，计算对流
            FOR_EACH_FACE{
                //半拉格朗日，倒退上一个时间的点
                glm::vec3 semi_LG;

                //防止越界
                if(i>-1 && i<grid.dim[0] && j>-1 && j<grid.dim[1] && k>-1 && k<grid.dim[2]){
                    //后边界，存储速度的x分量
                    semi_LG = grid.semiLagrangian(grid.getBack(i,j,k), Lagrangian2dPara::dt);
                    compute_X(i,j,k) = grid.getVelocityX(semi_LG);

                    //左边界，存储速度的y分量
                    semi_LG = grid.semiLagrangian(grid.getLeft(i,j,k), Lagrangian2dPara::dt);
                    compute_Y(i,j,k) = grid.getVelocityY(semi_LG);

                    //下边界，存储速度的z分量
                    semi_LG = grid.semiLagrangian(grid.getBottom(i,j,k), Lagrangian2dPara::dt);
                    compute_Z(i,j,k) = grid.getVelocityZ(semi_LG);
                }
                //中心点，存储密度和温度
                semi_LG = grid.semiLagrangian(grid.getCenter(i,j,k), Lagrangian2dPara::dt);
                compute_D(i,j,k) = grid.mD.interpolate(semi_LG);
                compute_T(i,j,k) = grid.mT.interpolate(semi_LG);

            }

            grid.mU=compute_X;
            grid.mV=compute_Y;
            grid.mW=compute_Z;
            grid.mD=compute_D;
            grid.mT=compute_T;

        }
		void Solver::compute_external_forces(MACGrid3d& grid){
            FOR_EACH_FACE{
                grid.mW(i,j,k)+=grid.getBoussinesqForce(grid.getCenter(i,j,k))*Lagrangian3dPara::dt;
            }
        }

		void Solver::projection(MACGrid3d& grid){
            //初始化压力为0
            std::vector<std::vector<std::vector<double>>>pressure(Eulerian3dPara::theDim3d[MACGrid3d::X],std::vector<std::vector<double>>(Eulerian3dPara::theDim3d[MACGrid3d::Y],std::vector<double>(Eulerian3dPara::theDim3d[MACGrid3d::Z],0.0)));

            double tar = 1e-7;//收敛值
            
            double rou = 1.225;//烟雾气体密度，假设烟雾的密度为常数，并接近于空气的密度
            double coeff = -(rou*grid.cellSize*grid.cellSize)/Lagrangian3dPara::dt;//基本的迭代系数

            for(int index=0; index<66; index++){
                auto compute_P = pressure;
                double maxR=0.0;

                FOR_EACH_CELL{
                    //忽略固体
                    if(grid.isSolidCell(i,j,k)){
                        continue;
                    }
                    double sum=0.0;
                    int count=0;

                    //对临近的6个cell计算
                    if(!grid.isSolidCell(i+1,j,k) && i+1<Eulerian3dPara::theDim3d[MACGrid3d::X]){
                        count++;
                        sum+=pressure[i+1][j][k];
                    }
                    if(!grid.isSolidCell(i,j+1,k) && j+1<Eulerian3dPara::theDim3d[MACGrid3d::Y]){
                        count++;
                        sum+=pressure[i][j+1][k];
                    }
                    if(!grid.isSolidCell(i,j,k+1) && k+1<Eulerian3dPara::theDim3d[MACGrid3d::Z]){
                        count++;
                        sum+=pressure[i][j][k+1];
                    }
                    if(!grid.isSolidCell(i-1,j,k) && i>0){
                        count++;
                        sum+=pressure[i-1][j][k];
                    }
                    if(!grid.isSolidCell(i,j-1,k) && j>0){
                        count++;
                        sum+=pressure[i][j-1][k];
                    }
                    if(!grid.isSolidCell(i,j,k-1) && k>0){
                        count++;
                        sum+=pressure[i][j][k-1];
                    }

                    if(count){
                        compute_P[i][j][k] = (sum+(grid.checkDivergence(i,j,k))*coeff)/count;
                    }

                    double R = compute_P[i][j][k]-pressure[i][j][k];
                    if(R<0)R=-R;
                    if(maxR<R){
                        maxR=R;
                    }

                }

                pressure=compute_P;
                if(maxR<tar){
                    break;//收敛
                }
            }

            FOR_EACH_CELL{
                if(!grid.isSolidCell(i,j,k)){
                    if(!grid.isSolidCell(i-1,j,k) && i>0){
                        double compute = pressure[i][j][k]-pressure[i-1][j][k];
                        grid.mU(i,j,k)-=(Lagrangian3dPara::dt/(grid.cellSize*rou))*compute;
                    }
                    if(!grid.isSolidCell(i,j-1,k) && j>0){
                        double compute = pressure[i][j][k]-pressure[i][j-1][k];
                        grid.mV(i,j,k)-=(Lagrangian3dPara::dt/(grid.cellSize*rou))*compute;
                    }
                    if(!grid.isSolidCell(i,j,k-1) && k>0){
                        double compute = pressure[i][j][k]-pressure[i][j][k-1];
                        grid.mW(i,j,k)-=(Lagrangian3dPara::dt/(grid.cellSize*rou))*compute;
                    }
                }
            }
        }

        void Solver::solve()
        {
            // TODO
            // Solves the fluid simulation by performing some steps, which may include:
            // 1. advection
            // 2. compute external forces
            // 3. projection
            // ...
            advection(mGrid);
            compute_external_forces(mGrid);
            projection(mGrid);
            
        }
    }
}

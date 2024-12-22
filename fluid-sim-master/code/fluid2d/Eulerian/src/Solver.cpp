#include "Eulerian/include/Solver.h"
#include "Configure.h"

namespace FluidSimulation
{
    namespace Eulerian2d
    {
        Solver::Solver(MACGrid2d &grid) : mGrid(grid)
        {
            mGrid.reset();
        }

        void Solver::advection(MACGrid2d& grid){
            auto compute_X = grid.mU;
            auto compute_Y = grid.mV;
            auto compute_D = grid.mD;
            auto compute_T = grid.mT;

            //遍历网格，计算对流
            FOR_EACH_LINE{
                //半拉格朗日，倒退上一个时间的点
                glm::vec2 semi_LG;

                //防止越界
                if(i>-1 && i<grid.dim[0] && j>-1 && j<grid.dim[1]){
                    //左边界，存储速度的x分量
                    semi_LG = grid.semiLagrangian(grid.getLeft(i,j), Lagrangian2dPara::dt);
                    compute_X(i,j) = grid.getVelocityX(semi_LG);

                    //下边界，存储速度的y分量
                    semi_LG = grid.semiLagrangian(grid.getBottom(i,j), Lagrangian2dPara::dt);
                    compute_Y(i,j) = grid.getVelocityY(semi_LG);
                }
                //中心点，存储密度和温度
                semi_LG = grid.semiLagrangian(grid.getCenter(i,j), Lagrangian2dPara::dt);
                compute_D(i,j) = grid.mD.interpolate(semi_LG);
                compute_T(i,j) = grid.mT.interpolate(semi_LG);

            }

            grid.mU=compute_X;
            grid.mV=compute_Y;
            grid.mD=compute_D;
            grid.mT=compute_T;

        }
		void Solver::compute_external_forces(MACGrid2d& grid){
            FOR_EACH_LINE{
                grid.mV(i,j)+=grid.getBoussinesqForce(grid.getCenter(i,j))*Lagrangian2dPara::dt;
            }
        }

		void Solver::projection(MACGrid2d& grid){
            //初始化压力为0
            std::vector<std::vector<double>>pressure(Eulerian2dPara::theDim2d[MACGrid2d::X],std::vector<double>(Eulerian2dPara::theDim2d[MACGrid2d::Y],0.0));

            double tar = 1e-7;//收敛值
            
            double rou = 1.3;//气体密度
            double coeff = -(rou*mGrid.cellSize*mGrid.cellSize)/Lagrangian2dPara::dt;//基本的迭代系数

            for(int index=0; index<66; index++){
                auto compute_P = pressure;
                double maxR=0.0;

                FOR_EACH_CELL{
                    //忽略固体
                    if(grid.isSolidCell(i,j)){
                        continue;
                    }
                    double sum=0.0;
                    int count=0;

                    //对临近的四个cell计算
                    if(!grid.isSolidCell(i+1,j) && i+1<Eulerian2dPara::theDim2d[MACGrid2d::X]){
                        count++;
                        sum+=pressure[i+1][j];
                    }
                    if(!grid.isSolidCell(i,j+1) && j+1<Eulerian2dPara::theDim2d[MACGrid2d::Y]){
                        count++;
                        sum+=pressure[i][j+1];
                    }
                    if(!grid.isSolidCell(i-1,j) && i>0){
                        count++;
                        sum+=pressure[i-1][j];
                    }
                    if(!grid.isSolidCell(i,j-1) && j>0){
                        count++;
                        sum+=pressure[i][j-1];
                    }

                    if(count){
                        compute_P[i][j] = (sum+(grid.checkDivergence(i,j))*coeff)/count;
                    }

                    double R = compute_P[i][j]-pressure[i][j];
                    if(maxR>R){
                        maxR=R;
                    }
                    else if(maxR>(-R)){
                        maxR=(-R);
                    }

                }

                pressure=compute_P;
                if(maxR<tar){
                    break;//收敛
                }
            }

            FOR_EACH_CELL{
                if(!mGrid.isSolidCell(i,j)){
                    if(!mGrid.isSolidCell(i-1,j) && i>0){
                        double compute = pressure[i][j]-pressure[i-1][j];
                        grid.mU(i,j)-=(Lagrangian2dPara::dt/(mGrid.cellSize*rou))*compute;
                    }
                    if(!mGrid.isSolidCell(i,j-1) && j>0){
                        double compute = pressure[i][j]-pressure[i][j-1];
                        grid.mV(i,j)-=(Lagrangian2dPara::dt/(mGrid.cellSize*rou))*compute;
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

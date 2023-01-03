//测试，调试，以及输出文件
#ifndef _diagnose_h
#define _diagnose_h
#include<string>
#include"poisson_solver.h"
#include"particles.h"

using namespace std;

void OutputData(string filename, vector<double> a);
vector<double> GetParticlesVX(Particles testp);
vector<double> GetParticlesX(Particles testp);

void OutputElectricPotentialGrids(string filename, PoissonSolverPeriodicBC_FFTW efg);
void OutputElectricFieldGrids(string filename, PoissonSolverPeriodicBC_FFTW efg);

#endif

#ifndef _poisson_solver_h
#define _poisson_solver_h

#include<eigen3/Eigen/Core>
#include <fftw3.h>
#include "particles.h"

using namespace std;
using namespace Eigen;

class PoissonSolverPeriodicBC_FFTW
{
  public:
    const int nx_grids;
    const double dx;
    const double L;
    VectorXd Ex;
    VectorXd phi;

    PoissonSolverPeriodicBC_FFTW(int _grids_num, double _grid_width);
    void Solve(VectorXd _charge);
    double GetEx(int x_idx);
    double GetPhi(int x_idx);
};

#endif

#ifndef _esfield_h
#define _esfield_h

#include <cmath>
#include<eigen3/Eigen/Core>
#include <map>
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include "particles.h"
#include "input.h"

using namespace std;
using namespace Eigen;

double LineInterpolation(double k, int kGrid, double gridWidth); //charge in this grid in k-th dimension

double QuadInterpolation(double x, int kGrid, double gridWidth, int index);

map<int, double> PartitionToGrid(double grid_width, int grids_num, double x, unsigned int interpolation);

// 3d array transform 1d
//int IndexFrom2dTo1d(int x, int y, int N); //N^3 is length of 1d array
//int IndexFrom1dTo2d(int i, int N);

class ElectricField
{
  public:
    const int nx_grids;
    const double dx;
    const double L;
    VectorXd Ex;
    VectorXd phi;

    ElectricField(int _grids_num, double _grid_width);
    void SolveEFieldOnGrids(VectorXd _charge);
    double GetEx(int x_idx);
    double GetPhi(int x_idx);
};

#endif

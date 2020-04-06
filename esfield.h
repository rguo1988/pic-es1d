#ifndef _esfield_h
#define _esfield_h

#include <cmath>
#include <map>
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include "particles.h"
#include "input.h"

using namespace std;

double LineInterpolation(double k, int kGrid, double gridWidth); //charge in this grid in k-th dimension

double QuadInterpolation(double x, int kGrid, double gridWidth, int index);
class GridPoint
{
  public:
    int x;
    int y;
    GridPoint(int xx = 0, int yy = 0);
};

bool operator<(const GridPoint& lhs, const GridPoint& rhs);

map<GridPoint, double> PartitionToGrid(double grid_width, int grids_num, PhaseSpace temp_ps, double x_min, unsigned int interpolation);

// 3d array transform 1d
int IndexFrom2dTo1d(int x, int y, int N); //N^3 is length of 1d array
GridPoint IndexFrom1dTo2d(int i, int N);

class ElectricFieldGrids: public UniversalParameters
{
  public:
    map<GridPoint, double> charge_dens_grid;
    map<GridPoint, double> elec_pot_grid;
    map<GridPoint, double> Ex_grid;

    //ElectricFieldGrids();
    void SetupEFieldOnGrids();
    double AppliedElecFields(double x, double t, double amplitude, double k, double w);
};

#endif

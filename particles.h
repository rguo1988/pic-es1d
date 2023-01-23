#ifndef _particles_h
#define _particles_h
#include<vector>
#include<string>
using namespace std;

class Particles
{
  public:
    const int num;//alpha-particle numbers
    const double q;//alpha-particle charge
    const double m;//alpha-particle mass
    const string name;//alpha-particle name
    vector<double> x;
    vector<double> vx;

    Particles(double N, double Q, double M, string name);
    void InitializeXV_Random(double (*Distribution)(double, double), double v_max, double L);
    void InitializeXV_Quiet(double (*DistributionX)(double), double (*DistributionV)(double), double L, double dx, int nx_grids);
};

#endif

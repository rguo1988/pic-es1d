#ifndef _particles_h
#define _particles_h
#include<vector>
#include<string>
using namespace std;

class PhaseSpace
{
  public:
    double x;
    double vx;

    PhaseSpace(double X = 0, double VX = 0);
    ~PhaseSpace(void);
};
class Particles
{
  public:
    const int num;//alpha-particle numbers
    const double q;//alpha-particle charge
    const double m;//alpha-particle mass
    const string name;//alpha-particle name
    vector<PhaseSpace> rv;

    Particles(double N, double Q, double M, string name);
};

#endif

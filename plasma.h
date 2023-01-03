#ifndef _plasma_h
#define _plasma_h

#include<gsl/gsl_rng.h>
#include<sys/timeb.h>
#include<eigen3/Eigen/Core>

#include"particles.h"
#include"poisson_solver.h"
#include"input.h"

using namespace std;
using namespace Eigen;

class PlasmaSystem: public Input
{
  public:
    PoissonSolverPeriodicBC_FFTW E;
    VectorXd charge;

    vector<double> Ek;
    vector<double> Ep;
    vector<double> Et;
    vector<double> T;
    vector<double> num_in_grid;

    PlasmaSystem();

    //print important parameters in simulation
    void PrintParameters() const;

    //trace energy during simulation
    void CalculateE();
    //calculate the temperature in every cell
    void CalculateT();

    //Setup Charge
    void SetupSpeciesChargeOnGrids();
    void SetupBackgroundChargeOnGrids();

    //BorisPusher
    PhaseSpace LangevinPusher(PhaseSpace last_rv, double gamma, double D, gsl_rng *r);
    void PushOneStep(int if_init);
    void Run();

    //excite another wave
    void ExicteWave(double A, double k);
};
#endif

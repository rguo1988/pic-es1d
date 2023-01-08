#ifndef _plasma_h
#define _plasma_h

#include"poisson_solver.h"
#include"input.h"

using namespace std;
using namespace Eigen;

class PlasmaSystem: public Input
{
  public:
    PoissonSolverPeriodicBC_FFTW E;
    VectorXd charge;
    VectorXd charge_background;

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
    void CalcBackgroundChargeOnGrids();

    //PhaseSpace LangevinPusher(PhaseSpace last_rv, double gamma, double D, gsl_rng *r);
    void Run();

    //excite another wave
    void ExicteWave(double A, double k);
};
#endif

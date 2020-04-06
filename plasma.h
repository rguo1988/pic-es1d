//等离子体系统
#ifndef _plasma_h
#define _plasma_h

#include<gsl/gsl_rng.h>
#include<sys/timeb.h>

#include"particles.h"
#include"esfield.h"
#include"bfield.h"
#include"input.h"

using namespace std;

class PlasmaSystem: public Input
{
  public:
    ElectricFieldGrids E;
    ConstMagneticField B;

    vector<double> Ek;
    vector<double> Ep;
    vector<double> Et;

    PlasmaSystem();

    //print important parameters in simulation
    void PrintParameters() const;

    //trace energy during simulation
    void CalculateE();
    
    //Setup Charge
    void ClearChargeOnGrids();
    void SetupSpeciesChargeOnGrids();
    void SetupBackgroundChargeOnGrids();

    //BorisPusher
    PhaseSpace BorisInitPusher(PhaseSpace init_rv, double fx, double mass);
    PhaseSpace BorisFinalPusher(PhaseSpace init_rv, double fx, double mass);
    PhaseSpace BorisPusher(PhaseSpace last_rv, double fx, double mass);
    void PushOneStep(int if_init);
    void Run();

    //excite another wave
    void ExicteWave(double A, double k);
};
#endif

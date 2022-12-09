#include "diagnose.h"
#include <iostream>
#include <iomanip>

using namespace std;

void OutputData(string filename, vector<double> a)
{
    ofstream ofile;
    ofile.open(filename.c_str());
    for (auto i : a)
    {
        ofile << setprecision(13) << i << endl;
    }
    ofile.close();
}
void OutputElectricPotentialGrids(string filename, PoissonSolverPeriodicBC_FFTW efg)
{
    ofstream ofile;
    ofile.open(filename.c_str());

    for(int i = 0; i < efg.nx_grids; i++)
    {
        ofile << efg.GetPhi(i) << endl;
    }
    ofile.close();
}
void OutputElectricFieldGrids(string filename, PoissonSolverPeriodicBC_FFTW efg)
{
    ofstream ofile;
    ofile.open(filename.c_str());

    for(int i = 0; i < efg.nx_grids; i++)
    {
        ofile << efg.GetEx(i)<< endl;
    }
    ofile.close();
}
vector<double> GetParticlesVX(Particles testp)
{
    vector<double> test;
    test.clear();
    for(auto i : testp.rv)
    {
        test.push_back(i.vx);
    }
    return test;
}
vector<double> GetParticlesX(Particles testp)
{
    vector<double> test;
    test.clear();
    for(auto i : testp.rv)
    {
        test.push_back(i.x);
    }
    return test;
}

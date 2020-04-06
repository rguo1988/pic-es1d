#include "diagnose.h"
#include <iostream>
#include <iomanip>

using namespace std;

double OutputTotCharge(ElectricFieldGrids cdg)
{
    double tempq = 0.0;
    for(auto i : cdg.charge_dens_grid)
    {
        tempq += i.second;
    }
    return tempq;
}
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
void OutputNetCharge(ElectricFieldGrids efg)
{
    double tot_charge = 0.0;
    //for(int i = 0; i < efg.grids_num; i++)
    for(auto rho : efg.charge_dens_grid)
    {
        tot_charge += rho.second;
    }
    if(tot_charge == 0)
        cout << "Net Charge is zero!" << endl;
    else
        cout << "Net Charge is NOT zero!!!" << tot_charge << endl;
}
void OutputChargeGrids(string filename, ElectricFieldGrids cdg)
{
    ofstream ofile;
    ofile.open(filename.c_str());

    for(auto g : cdg.charge_dens_grid)
    {
        ofile << g.first.x << "  " << g.second << endl;
    }
    ofile.close();
}
void OutputElectricPotentialGrids(string filename, ElectricFieldGrids efg)
{
    ofstream ofile;
    ofile.open(filename.c_str());

    for(auto g : efg.elec_pot_grid)
    {
        ofile << g.first.x << "  " << g.second << endl;
    }
    ofile.close();
}
void OutputElectricFieldGrids(string filename, ElectricFieldGrids efg)
{
    ofstream ofile;
    ofile.open(filename.c_str());

    for(auto g : efg.Ex_grid)
    {
        GridPoint tgp(g.first.x);
        ofile << g.first.x << "  " << g.second << endl;
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

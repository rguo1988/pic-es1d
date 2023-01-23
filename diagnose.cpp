#include<iomanip>
#include<fstream>
#include"diagnose.h"
#include"plasma.h"

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

void PlasmaSystem::DiagnoseDistribution(string p)
{

    for(auto particles_a : species)
    {
        string filenameVX = data_path + particles_a.name + "_vx_data" + p;
        string filenameX = data_path + particles_a.name + "_x_data" + p;
        //output particles x v
        OutputData(filenameVX, particles_a.vx);
        OutputData(filenameX,  particles_a.x);
    }
}

void PlasmaSystem::DiagnoseEnergy()
{
    double tempEk = 0.0;
    int particles_tot_num = 0;
    for (auto particles_a : species)
    {
        for(int i = 0; i < particles_a.num; i++)
        {
            tempEk += 0.5 * particles_a.m * pow(particles_a.vx[i], 2);
        }
        particles_tot_num += particles_a.num;
    }
    tempEk /= particles_tot_num;
    Ek.push_back(tempEk);

    double tempEp = 0.0;
    for(int i = 0; i < nx_grids; i++)
    {
        tempEp += 0.5 * pow(E.GetE(i), 2) * dx;
    }
    tempEp /= particles_tot_num;
    Ep.push_back(tempEp);
    Et.push_back(tempEk + tempEp);
}

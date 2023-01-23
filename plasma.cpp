//class plasma realization
#include"plasma.h"
#include"partition1d.h"
#include"input.h"
#include"diagnose.h"
#include<iostream>
#include<iomanip>
#include<fstream>

using namespace std;

void PlasmaSystem::DiagnoseTemperature()
{
    T.clear();
    T.resize(nx_grids, 0.0);
    num_in_grid.clear();
    num_in_grid.resize(nx_grids, 0.0);

    for (auto particles_a : species)
    {
        for(int i = 0; i < particles_a.num; i++)
        {
            int idx_grid = floor(particles_a.x[i] / dx);
            if (idx_grid == nx_grids)
                idx_grid = 0;
            T[idx_grid] += particles_a.m * pow(particles_a.vx[i], 2);
            num_in_grid[idx_grid]++;
        }
    }
    for(int i = 0; i < nx_grids; i++)
    {
        T[i] /= num_in_grid[i];
    }
}
PlasmaSystem::PlasmaSystem():
    E(nx, dx)
{
    Ek.clear();
    Ep.clear();
    Et.clear();
    charge.resize(nx_grids);
    charge_background.resize(nx_grids);
}

//PhaseSpace PlasmaSystem::LangevinPusher(PhaseSpace last_rv, double gamma, double D, gsl_rng* r)
//{
//int idx = floor(last_rv.x / dx);
//if (idx == nx_grids)
//idx = 0;
//D = gamma * T[idx] / m_i;

//double tempx = last_rv.x + last_rv.vx * dt;
//double tempv = last_rv.vx - gamma * last_rv.vx * dt + sqrt(D * dt ) * gsl_ran_gaussian(r, sqrt(2.0));
//PhaseSpace next_rv(tempx, tempv);
//return next_rv;
//return 0.0;
//}

void PlasmaSystem::ExicteWave(double A, double k)
{
    for(int i = 0; i < species[0].num; i++)
    {
        double xplus = A * cos(k * species[0].x[i]);
        species[0].x[i] += xplus;
        while(species[0].x[i] < 0.0)
            species[0].x[i] += L;
        while(species[0].x[i] >= L)
            species[0].x[i] -= L;
    }
    cout << "Wave is Excited!" << endl
         << "A= " << A
         << "; k= " << k
         << "; k*l_D=" << k*lambda_D << endl;
}

void PlasmaSystem::Run()
{
    Initialize();
    PrintParameters();
    CalcBackgroundChargeOnGrids();

    //main loop
    for(int n = 0; n < maxsteps + 1; n++)
    {
        DiagnoseTemperature();
        int percent = 100 * n / (maxsteps);
        //print running process
        if(percent % 5 == 0)
        {
            cout << "\r" << "  Calculation Process:" << percent << "%" << flush;
        }

        //diagnose plasma every timestep
        if(n % data_steps == 0)
        {
            string p = to_string(n / data_steps);
            DiagnoseDistribution(p);
        }

        //calculate E
        charge.setZero();
        SetupSpeciesChargeOnGrids();
        //add neutral background charge density
        charge += charge_background;
        E.Solve(charge);

        //PushOneStep;
        for(auto &particles_a : species)
        {
            #pragma omp parallel for
            for(int j = 0; j < particles_a.num; j++)
            {
                map<int, double> partition_contrib = PartitionToGrid(dx, nx_grids, particles_a.x[j], 2); //linear interpolation //ec scheme
                double fex = 0.0;
                for(auto g : partition_contrib)
                {
                    fex += particles_a.q * E.GetE(g.first) * g.second * dx;
                }
                if(n == 0)
                {
                    particles_a.vx[j] += 0.5 * fex * dt / particles_a.m;
                    particles_a.x[j]  += particles_a.vx[j] * dt;
                }
                else
                {
                    particles_a.vx[j] += fex * dt / particles_a.m;
                    particles_a.x[j]  += particles_a.vx[j] * dt;
                }
                //period condition
                while(particles_a.x[j] < 0.0)
                {
                    particles_a.x[j] += L;
                }
                while(particles_a.x[j] >= L)
                {
                    particles_a.x[j] -= L;
                }
            }
        }
        DiagnoseEnergy();
    }
    //output energy evolution
    OutputData(data_path + "/tot_energy", Et);
    OutputData(data_path + "/kin_energy", Ek);
    OutputData(data_path + "/pot_energy", Ep);

    cout << endl << "  Complete!" << endl;
}

void PlasmaSystem::PrintParameters()
{
    cout << "----------------------------------------------------------------------" << endl;
    cout << "  1D ES PIC Simulation Start!" << endl;
    cout << "----------------------------------------------------------------------" << endl;
    cout << "  Simulation Parameters:" << endl;
    cout << "     L = " << left << setw(7) << setprecision(4)  << L
         << "     k = " << left << setw(7) << k
         << "    nx = " << left << setw(7) << nx
         << "    dx = " << left << setw(7) << setprecision(4) << dx << endl;

    cout << "   w_p = " << left << setw(7) << setprecision(4)  << w_pe
         << "   l_D = " << left << setw(7) << lambda_D << endl;

    cout << " Steps = " << left << setw(7) << maxsteps
         << "    dt = " << left << setw(7)  << dt
         << "  Time = " << left << setw(7)  << maxsteps*dt << endl;
    cout << "----------------------------------------------------------------------" << endl;
    if(dx > lambda_D)
    {
        cout << "WARNING: DebyeL is NOT satisfied!" << endl;
        cout << "----------------------------------------------------------------------" << endl;
    }
    for(auto particles_a : species)
    {
        cout << "  Species: " << particles_a.name << endl;
        cout << "     m = " << left << setw(7) << particles_a.m
             << "N/Cell = " << left << setw(7) << particles_a.num / nx_grids
             << "     N = " << left << setw(7) << particles_a.num
             << "     q = " << left << setw(7) << setprecision(4) << particles_a.q << endl;
    }
    cout << "----------------------------------------------------------------------" << endl;

    // print special information in this experiment
    PrintSpecialInformation();

    cout << "  Data: " << endl;
    if(!if_continued)
        cout << "  Evolving from a NEW initial state!" << endl;
    else
        cout << "  Evolving from a CONTINUED state!" << endl;
    cout << "  data_num = " << data_num << endl;
    cout << "----------------------------------------------------------------------" << endl;
}

void PlasmaSystem::CalcBackgroundChargeOnGrids()
{
    double net_charge = 0.0;
    for(auto p : species)
    {
        net_charge += p.num * p.q;
    }
    double normalization = 0.0;
    VectorXd temp_rho(nx_grids);
    for(int i = 0; i < nx_grids; i++)
    {
        double x = i * dx;
        temp_rho[i] = GetBackgroundDensity(x);
        normalization += temp_rho[i];
    }
    charge_background = net_charge * temp_rho / dx / normalization;
}

void PlasmaSystem::SetupSpeciesChargeOnGrids()
{
    //CIC 2-order interpolation
    for(auto p : species)
    {
        for(int i = 0; i < p.num; i++)
        {
            map<int, double> partition_contrib = PartitionToGrid(dx, nx_grids, p.x[i], 2); //interpolation
            for(auto g : partition_contrib)
            {
                charge[g.first] += p.q * g.second;// / pow(gridWidth, 1);
            }
        }
    }
}

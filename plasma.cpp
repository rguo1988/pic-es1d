//class plasma realization
#include"plasma.h"
#include"partition1d.h"
#include"input.h"
#include"diagnose.h"
#include<iostream>
#include<iomanip>
#include<fstream>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

using namespace std;

void PlasmaSystem::CalculateE()
{
    double tempEk = 0.0;
    int particles_tot_num = 0;
    for (auto particles_a : species)
    {
        for(auto prv : particles_a.rv)
        {
            tempEk += 0.5 * particles_a.m * (prv.vx * prv.vx);
        }
        particles_tot_num += particles_a.num;
    }
    tempEk /= particles_tot_num;
    Ek.push_back(tempEk);

    double tempEp = 0.0;
    for(int i = 0; i < nx_grids; i++)
    {
        tempEp += 0.5 * pow(E.GetEx(i), 2) * dx;
    }
    tempEp /= particles_tot_num;
    Ep.push_back(tempEp);
    Et.push_back(tempEk + tempEp);
}
void PlasmaSystem::CalculateT()
{
    T.clear();
    T.resize(nx_grids, 0.0);
    num_in_grid.clear();
    num_in_grid.resize(nx_grids, 0.0);

    for (auto particles_a : species)
    {
        for(auto prv : particles_a.rv)
        {
            int idx_grid = floor(prv.x / dx);
            if (idx_grid == nx_grids)
                idx_grid = 0;
            T[idx_grid] += particles_a.m * (prv.vx * prv.vx);
            num_in_grid[idx_grid]++;
        }
    }
    for(int i = 0; i < nx_grids; i++)
    {
        T[i] /= num_in_grid[i];
    }
}
PlasmaSystem::PlasmaSystem():
    B(0, 0, 0), E(nx_grids, dx)
{
    Ek.clear();
    Ep.clear();
    Et.clear();
    charge.resize(nx_grids);
}

PhaseSpace PlasmaSystem::LangevinPusher(PhaseSpace last_rv, double gamma, double D, gsl_rng* r)
{
    //int idx = floor(last_rv.x / dx);
    //if (idx == nx_grids)
    //idx = 0;
    //D = gamma * T[idx] / m_i;

    //double tempx = last_rv.x + last_rv.vx * dt;
    //double tempv = last_rv.vx - gamma * last_rv.vx * dt + sqrt(D * dt ) * gsl_ran_gaussian(r, sqrt(2.0));
    //PhaseSpace next_rv(tempx, tempv);
    //return next_rv;
    return 0.0;
}
void PlasmaSystem::PushOneStep(int if_init)
{
    /*
    // generating rnd number
    struct timeb time_seed;
    ftime(&time_seed);
    gsl_rng_default_seed = (time_seed.time * 1000 + time_seed.millitm);
    gsl_rng *r;
    r = gsl_rng_alloc(gsl_rng_default);

    for(auto &particles_a : species)
    {
        #pragma omp parallel for
        for(int j = 0; j < particles_a.num; j++)
        {
            map<int, double> partition_contrib = PartitionToGrid(dx, nx_grids, particles_a.rv[j].x, 2); //linear interpolation //ec scheme
            double fex = 0.0;
            for(auto g : partition_contrib)
            {
                fex += particles_a.q * E.GetEx(g.first) * g.second * dx;
            }
            if(if_init == 0)
            {
                particles_a.rv[j].vx += 0.5 * fex * dt / particles_a.m;
                particles_a.rv[j].x += particles_a.rv[j].vx * dt;
            }
            else
            {
                particles_a.rv[j].vx += fex * dt / particles_a.m;
                particles_a.rv[j].x += particles_a.rv[j].vx * dt;
                //particles_a.rv[j] = LangevinPusher(particles_a.rv[j], gamma, D, r);
            }
            //period condition
            while(particles_a.rv[j].x < 0.0)
            {
                particles_a.rv[j].x += L;
            }
            while(particles_a.rv[j].x >= L)
            {
                particles_a.rv[j].x -= L;
            }
        }
    }
    */
}

void PlasmaSystem::ExicteWave(double A, double k)
{
    for(int i = 0; i < species[0].num; i++)
    {
        double xplus = A * cos(k * species[0].rv[i].x);
        species[0].rv[i].x += xplus;
        while(species[0].rv[i].x < 0.0)
            species[0].rv[i].x += L;
        while(species[0].rv[i].x >= L)
            species[0].rv[i].x -= L;
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
    PrintSpecialInformation();

    //main loop
    for(int n = 0; n < maxsteps + 1; n++)
    {
        CalculateT();
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

            for(auto particles_a : species)
            {
                string filenameV = data_path + particles_a.name + "v_data" + p;
                string filenameX = data_path + particles_a.name + "x_data" + p;
                //output particles x v
                OutputData(filenameV, GetParticlesVX(particles_a));
                OutputData(filenameX,  GetParticlesX(particles_a));
            }
        }

        //calculate E
        charge.resize(nx_grids);
        charge.setZero();
        SetupSpeciesChargeOnGrids();
        SetupBackgroundChargeOnGrids();
        E.Solve(charge);

        //PushOneStep;
        //struct timeb time_seed;
        //ftime(&time_seed);
        //gsl_rng_default_seed = (time_seed.time * 1000 + time_seed.millitm);
        //gsl_rng *r;
        //r = gsl_rng_alloc(gsl_rng_default);

        for(auto &particles_a : species)
        {
            #pragma omp parallel for
            for(int j = 0; j < particles_a.num; j++)
            {
                map<int, double> partition_contrib = PartitionToGrid(dx, nx_grids, particles_a.rv[j].x, 2); //linear interpolation //ec scheme
                double fex = 0.0;
                for(auto g : partition_contrib)
                {
                    fex += particles_a.q * E.GetEx(g.first) * g.second * dx;
                }
                if(n == 0)
                {
                    particles_a.rv[j].vx += 0.5 * fex * dt / particles_a.m;
                    particles_a.rv[j].x += particles_a.rv[j].vx * dt;
                }
                else
                {
                    particles_a.rv[j].vx += fex * dt / particles_a.m;
                    particles_a.rv[j].x += particles_a.rv[j].vx * dt;
                    //particles_a.rv[j] = LangevinPusher(particles_a.rv[j], gamma, D, r);
                }
                //period condition
                while(particles_a.rv[j].x < 0.0)
                {
                    particles_a.rv[j].x += L;
                }
                while(particles_a.rv[j].x >= L)
                {
                    particles_a.rv[j].x -= L;
                }
            }
        }
        CalculateE();
    }
    //output energy evolution
    OutputData(data_path + "/tot_energy", Et);
    OutputData(data_path + "/kin_energy", Ek);
    OutputData(data_path + "/pot_energy", Ep);

    cout << endl << "  Complete!" << endl;
}

void PlasmaSystem::PrintParameters() const
{
    cout << "--------------------------------------------" << endl;
    cout << "  PIC Simulation Start!" << endl;
    cout << "--------------------------------------------" << endl;
    cout << "  Simulation Parameters:" << endl;
    cout << setw(13) << "  Length"
         << setw(13) << "       k"
         << setw(13) << "nx_grids"
         << setw(13) << "Lambda_D"
         << setw(13) << setprecision(6) << "      dx" << endl;

    cout << setw(13) << L
         << setw(13) << k
         << setw(13) << nx_grids
         << setw(13) << lambda_D
         << setw(13) << dx << endl;

    cout << setw(13) << "MaxSteps"
         << setw(13) << "      dt"
         << setw(13) << "    Time" << endl;

    cout << setw(13) << maxsteps
         << setw(13) << dt
         << setw(13) << maxsteps*dt << endl;

    cout << "--------------------------------------------" << endl;
    if(dx > lambda_D)
    {
        cout << "WARNING: DebyeL is NOT satisfied!" << endl;
        cout << "--------------------------------------------" << endl;
    }
    for(auto particles_a : species)
    {
        cout << "  Species: " << particles_a.name << endl;
        cout  << setw(8) << "N = "  << particles_a.num
              << setw(15) << "N/Cell = " << particles_a.num / nx_grids
              << setw(8) << "q = " << setprecision(6) << particles_a.q
              << setw(8) << "m = " << particles_a.m << endl;
    }
    cout << "--------------------------------------------" << endl;
    cout << "  Data: " << endl;
    if(!if_continued)
        cout << "  Evolving from a NEW initial state!" << endl;
    else
        cout << "  Evolving from a CONTINUED state!" << endl;
    cout << "  data_num = " << data_num << endl;
    cout << "--------------------------------------------" << endl;

}

void PlasmaSystem::SetupBackgroundChargeOnGrids()
{
    double net_charge = 0.0;
    for(auto p : species)
    {
        net_charge += p.num * p.q;
    }
    double normalization = 0.0;
    vector<double> rho(nx_grids, 0.0);
    for(int i = 0; i < nx_grids; i++)
    {
        double x = i * dx;
        rho[i] = GetBackgroundDensity(x);
        normalization += rho[i];
    }
    for(int i = 0; i < nx_grids; i++)
    {
        rho[i] /= normalization;
        charge[i] -= net_charge * rho[i] / dx;
    }
}

void PlasmaSystem::SetupSpeciesChargeOnGrids()
{
    //CIC 2-order interpolation
    for(auto p : species)
    {
        for(int i = 0; i < p.num; i++)
        {
            map<int, double> partition_contrib = PartitionToGrid(dx, nx_grids, p.rv[i].x, 2); //interpolation
            for(auto g : partition_contrib)
            {
                charge[g.first] += p.q * g.second;// / pow(gridWidth, 1);
            }
        }
    }
}

//class plasma realization
#include"plasma.h"
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
            tempEk += particles_a.m * (prv.vx * prv.vx) / 2;
        }
        particles_tot_num += particles_a.num;
    }
    tempEk /= particles_tot_num;
    Ek.push_back(tempEk);

    double tempEp = 0.0;
    for(auto s : E.charge_dens_grid)
    {
        tempEp += E.elec_pot_grid[s.first] * s.second * pow(grid_width, 1) / 2;
    }
    tempEp /= particles_tot_num;
    Ep.push_back(tempEp);
    Et.push_back(tempEk + tempEp);
}
void PlasmaSystem::CalculateT()
{
    for (auto particles_a : species)
    {
        for(auto prv : particles_a.rv)
        {

            tempEk += particles_a.m * (prv.vx * prv.vx) / 2;
        }
        particles_tot_num += particles_a.num;
    }
}
PlasmaSystem::PlasmaSystem():
    B(0, 0, 0)
{

    Ek.clear();
    Ep.clear();
    Et.clear();
}

PhaseSpace PlasmaSystem::BorisInitPusher(PhaseSpace init_rv, double fx, double mass)
{
    double v_half = init_rv.vx + fx * dt / mass / 2;
    double x = init_rv.x + v_half * dt;
    PhaseSpace r_half_v(x, v_half);
    return r_half_v;
}

PhaseSpace PlasmaSystem::BorisFinalPusher(PhaseSpace init_rv, double fx, double mass)
{
    double v_half = init_rv.vx + fx * dt / mass / 2;
    PhaseSpace r_half_v(init_rv.x, v_half);
    return r_half_v;
}

PhaseSpace PlasmaSystem::BorisPusher(PhaseSpace last_rv, double fx, double mass)
{
    double vx_next = last_rv.vx + fx * dt / mass;
    double x_next = last_rv.x + dt * vx_next;

    PhaseSpace next_rv(x_next, vx_next);
    return next_rv;
}

PhaseSpace PlasmaSystem::LangevinPusher(PhaseSpace last_rv, double gamma, double D, gsl_rng* r)
{
    // dv= - gamma * v * dt + sqrt(D * dt) * w;
    double tempx = last_rv.x + last_rv.vx * dt;
    double tempv = last_rv.vx - gamma * last_rv.vx * dt + sqrt(D * dt ) * gsl_ran_gaussian(r, sqrt(2.0));
    PhaseSpace next_rv(tempx, tempv);
    return next_rv;
}
void PlasmaSystem::PushOneStep(int if_init)
{
    // generating rnd number
    struct timeb time_seed;
    ftime(&time_seed);
    gsl_rng_default_seed = (time_seed.time * 1000 + time_seed.millitm);
    gsl_rng *r;
    r = gsl_rng_alloc(gsl_rng_default);
    //double w = gsl_ran_gaussian(r, sqrt(2.0));

    for(auto &particles_a : species)
    {
        #pragma omp parallel for
        for(int j = 0; j < particles_a.num; j++)
        {
            map<GridPoint, double> partition_contrib = PartitionToGrid(grid_width, grids_num, particles_a.rv[j], x_min, 2); //linear interpolation //ec scheme
            double fex = 0.0;
            for(auto g : partition_contrib)
            {
                fex += particles_a.q * E.Ex_grid[g.first] * g.second * grid_width;
            }

            if(if_init == 0)
            {
                particles_a.rv[j] = BorisInitPusher(particles_a.rv[j], fex, particles_a.m);
            }
            else
            {
                particles_a.rv[j] = BorisPusher(particles_a.rv[j], fex, particles_a.m);
                //particles_a.rv[j] = LangevinPusher(particles_a.rv[j], gamma, D, r);
            }
            //period condition
            while(particles_a.rv[j].x < x_min)
            {
                particles_a.rv[j].x += L;
            }
            while(particles_a.rv[j].x >= x_max)
            {
                particles_a.rv[j].x -= L;
            }
        }
    }
}

void PlasmaSystem::ExicteWave(double A, double k)
{
    for(int i = 0; i < species[0].num; i++)
    {
        double xplus = A * cos(k * species[0].rv[i].x);
        //double vxplus = A * k * sin(k * electrons.rv[i].x);
        species[0].rv[i].x += xplus;
        //electrons.rv[i].vx += vxplus;
        while(species[0].rv[i].x < x_min)
            species[0].rv[i].x += L;
        while(species[0].rv[i].x >= x_max)
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
    for(int n = 0; n < maxsteps; n++)
    {
        int percent = 100 * n / (maxsteps - 1);
        //print running process
        if(percent % 5 == 0)
        {
            cout << "\r" << "Calculation Process:" << percent << "%" << flush;
        }

        //diagnose plasma every timestep
        if(n % data_steps == 0)
        {
            string p = to_string(n / data_steps);

            for(auto particles_a : species)
            {
                string filenameV = data_path + particles_a.name + "v_data" + p;
                string filenameX = data_path + particles_a.name + "x_data" + p;
                OutputData(filenameV, GetParticlesVX(particles_a));
                OutputData(filenameX,  GetParticlesX(particles_a));
            }
        }

        //calculate E
        ClearChargeOnGrids();
        SetupSpeciesChargeOnGrids();
        SetupBackgroundChargeOnGrids();
        E.SetupEFieldOnGrids();
        PushOneStep(n);
        CalculateE();
    }
    //output energy evolution
    OutputData(data_path + "/tot_energy", Et);
    OutputData(data_path + "/kin_energy", Ek);
    OutputData(data_path + "/pot_energy", Ep);

    cout << endl << "Complete!" << endl;
}

void PlasmaSystem::PrintParameters() const
{
    cout << "--------------------------------------------" << endl;
    cout << "PIC Simulation Start!" << endl;
    cout << "--------------------------------------------" << endl;
    cout << "Simulation Parameters:" << endl;
    cout << setw(10) << "Length"
         << setw(15) << "Grids Num"
         << setw(15) << "Grid Width"
         << setw(15) << "P Per Cell" << endl;

    cout << setw(10) << L
         << setw(15) << grids_num
         << setw(15) << grid_width
         << setw(15) << N_e / grids_num << endl;

    cout << setw(10) << "Max Step"
         << setw(15) << "dt"
         << setw(15) << "Time"
         << setw(15) << "Total Time" << endl;

    cout << setw(10) << maxsteps
         << setw(15) << dt
         << setw(15) << maxsteps*timestep_condition
         << setw(15) << maxsteps*timestep_condition + time_ran << endl;
    cout << "--------------------------------------------" << endl;
    if(grid_width > lambda_D)
    {
        cout << "WARNING: DebyeL is NOT satisfied!" << endl;
        cout << "--------------------------------------------" << endl;
    }
    for(auto particles_a : species)
    {
        cout << "Species: " << particles_a.name << endl;
        cout  << setw(5) << "N = "  << particles_a.num
              << setw(10) << "q = " << particles_a.q
              << setw(10) << "m = " << particles_a.m << endl;
    }
    cout << "--------------------------------------------" << endl;
    cout << "Data: " << endl;
    if(!if_continued)
        cout << "Evolving from a NEW initial state!" << endl;
    else
        cout << "Evolving from a CONTINUED state!" << endl;
    cout << "data_num = " << data_num << endl;
    cout << "--------------------------------------------" << endl;

}

void PlasmaSystem::ClearChargeOnGrids()
{
    E.charge_dens_grid.clear();
    for (int gNumX = 0; gNumX < grids_num; gNumX++)
    {
        GridPoint temp_gp(gNumX);
        E.charge_dens_grid.insert(make_pair(temp_gp, 0.0));
    }
}

void PlasmaSystem::SetupBackgroundChargeOnGrids()
{
    double net_charge = 0.0;
    for(auto p : species)
    {
        net_charge += p.num * p.q;
    }
    double normalization = 0.0;
    vector<double> rho(grids_num, 0.0);
    for(int i = 0; i < grids_num; i++)
    {
        double x = i * grid_width + x_min;
        rho[i] = GetBackgroundIonDensity(x);
        normalization += rho[i];
    }
    for(int i = 0; i < grids_num; i++)
    {
        rho[i] /= normalization;
        E.charge_dens_grid[i] -= net_charge * rho[i] / grid_width;
    }
}

void PlasmaSystem::SetupSpeciesChargeOnGrids()
{
    //CIC 2-order interpolation
    //#pragma omp parallel for
    for(auto p : species)
    {
        for(int i = 0; i < p.num; i++)
        {
            map<GridPoint, double> partition_contrib = PartitionToGrid(grid_width, grids_num, p.rv[i], x_min, 2); //interpolation
            for(auto g : partition_contrib)
            {
                E.charge_dens_grid[g.first] += p.q * g.second;// / pow(gridWidth, 1);
            }
        }
    }
}

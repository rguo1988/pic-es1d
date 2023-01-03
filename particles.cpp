//
#include "particles.h"
#include<gsl/gsl_rng.h>
#include<ctime>
#include<cmath>
#include<iostream>

using namespace std;

PhaseSpace::PhaseSpace(double X, double VX)
{
    x = X;
    vx = VX;
}
PhaseSpace::~PhaseSpace(void)
{}

Particles::Particles(double _n, double _q, double _m, string _name):
    num(_n), q(_q),  m(_m), name(_name)
{
    rv.resize(num, 0.0);
}
void Particles::InitializeXV_Random(double (*Distribution)(double, double), double v_max, double L)
{
    //set random number
    gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    gsl_rng *r;
    r = gsl_rng_alloc(gsl_rng_default);

    //initialize particles
    double max_probability_density = Distribution(0.0, 0.0);
    for(int i = 0; i < num; i++)
    {
        double temp_vx = gsl_rng_uniform(r) * 2 * v_max - v_max;
        double temp_x = gsl_rng_uniform(r) * L;

        while(gsl_rng_uniform(r) * max_probability_density > Distribution(temp_x, temp_vx))
        {
            temp_x = gsl_rng_uniform(r) * L;
            temp_vx = gsl_rng_uniform(r) * 2 * v_max - v_max;
        }
        rv[i].x = temp_x;
        rv[i].vx = temp_vx;
    }


}
void Particles::InitializeXV_Quiet(double (*DistributionX)(double), double (*DistributionV)(double), double L, double dx, int nx_grids)
{
    //quiet start
    int Nx[nx_grids];
    int temp_N = 0;//particle numbers that have been assigned
    double dv = 1e-7;
    bool minus_or_plus = 0;

    for(int i = 0; i < nx_grids; i++)
    {
        //print assigned process
        int process = 100 * i / (nx_grids - 1);
        if(process % 5 == 0)
            cout << "\r" << "  Assignning " << name << " (Quiet Start): " << process << "%" << flush;

        double temp_x = (i + 0.5) * dx;
        Nx[i] = round(1.0 * num / nx_grids * DistributionX(temp_x));
        for(int j = temp_N; j < temp_N + Nx[i]; j++)
            rv[j].x = temp_x;
        double u1 = 0.0;
        for(int k = temp_N; k < temp_N + Nx[i] / 2.0; k++)
        {
            double integral = 0.0;
            int l = 0;
            for(; l < 5 / dv; l++)
            {
                integral += Nx[i] * DistributionV(u1 + l * dv) * dv;
                if(integral > 1.0)
                    break;
            }
            double temp_v = u1 + 0.5 * l * dv;
            rv[k].vx = temp_v;
            int k_inverse = 2 * temp_N + Nx[i] - 1 - k;
            if(k != k_inverse)
                rv[k_inverse].vx = -temp_v;
            else
            {
                rv[k_inverse].vx = (minus_or_plus) ? -temp_v : temp_v;
                minus_or_plus = !minus_or_plus;
            }
            u1 += l * dv;
        }
        temp_N += Nx[i];
    }
    if(temp_N != num)
        cout << "Quiet start ERROR!" << endl;
    cout << endl;
}

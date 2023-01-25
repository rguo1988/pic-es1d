#ifndef _input_h
#define _input_h
#include<cmath>
#include"particles.h"
#include<string>
#include<iostream>
using namespace std;

class Input
{
  public:
    //title
    const string title = "Landau damping of Ion-acoustic waves (using Boltzmann electrons as background)";
    //if continue from data
    const bool if_continued = 0;

    //simulation box
    static constexpr double k = 1.0;
    static constexpr double L = 2.0 * M_PI / k;
    const double v_max = 5.0;
    const double vx_width = 2.0 * v_max;

    static const int nx = 512;
    static const int nx_grids = nx - 1;
    const double dx = L / nx_grids;

    //speices
    static constexpr double m_e = 1.0;
    static constexpr double NePerCell = 500;
    static constexpr double T_e = 1;
    const double N_e = NePerCell * nx_grids;
    const double n_e_aver = N_e / L;
    const double q_e = -sqrt(1.0 / n_e_aver);
    const double w_pe = sqrt(1.0 / m_e);//set nq^2=1
    const double lambda_De = sqrt(T_e);

    static constexpr double m_i = 100.0;
    static constexpr double NiPerCell = 500;
    static constexpr double T_i = 0.1;
    const double N_i = NiPerCell * nx_grids;
    const double n_i_aver = N_i / L;
    const double q_i = sqrt(1.0 / n_i_aver);
    const double w_pi = sqrt(1.0 / m_i);//set nq^2=1
    const double lambda_Di = sqrt(T_i);

    const double w_p = sqrt(w_pe * w_pe + w_pi*w_pi);
    const double lambda_D = 1.0 / sqrt(1.0 / T_e + 1.0 / T_i);

    //time parameters
    const int maxsteps = 1500;
    const int time_ran = 0;
    const double timestep_condition = 0.1;
    //const double dt = timestep_condition / w_p;
    const double dt = 0.1;

    //special settings
    static constexpr double d = 0.01 * L;

    //data path
    const string data_path = "./data/";
    const int data_steps = maxsteps;
    const int data_num = maxsteps / data_steps + 1;

    vector<Particles> species;

    static double GetIonInitDistrib(double x, double v)
    {
        double f = sqrt( m_i / (2 * M_PI * T_i) ) * exp(-0.5 * m_i * v * v / T_i);
        return f / L;
    }

    // use equilibrium electrons as background ne = n0*exp(e\phi/T_e) = n0 (1+e\phi/T_e)
    double GetBackgroundDensity(double x, double phi_i)
    {
        // return 1.0 + q_i * phi_i / T_e;
        return exp(q_i * phi_i / T_e);
    }

    void Initialize()
    {
        Particles ions(N_i, q_i, m_i, "ions");
        cout << "----------------------------------------------------------------------" << endl;
        ions.InitializeXV_Random(GetIonInitDistrib, v_max, L);
        for(int i = 0; i < ions.num; i++)
        {
            ions.x[i] += d * cos(k * ions.x[i]);
            //apply periodic boundary conditions
            if(ions.x[i] > L)
                ions.x[i] -= L;
            if(ions.x[i] < 0)
                ions.x[i] += L;
        }
        species.push_back(ions);
    }
    void PrintSpecialInformation()
    {
        cout << "  " << title << endl;
        cout << "  " << "d = " << d << endl;
        cout << "----------------------------------------------------------------------" << endl;
    }
};

#endif

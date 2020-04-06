#ifndef _input_h
#define _input_h
#include<string>
#include<cmath>
#include"particles.h"

class Input
{
  public:
    //if the simulation continued from last calculation
    const bool if_continued = 0;
    //electron mass
    static constexpr double m_e = 1.0;
    //electron number
    static constexpr double N_e = 1e5;
    //electron temperature
    static constexpr double T_e = 1.0;

    //ion charge
    const double q_i = 0;
    //ion mass
    const double m_i = 0;
    //ion number
    const double N_i = 0;
    //ion temperature
    const double T_i = 0;

    //configuration space from x_min to x_max, from vx_min to vx_max
    const double x_min = -10.0;
    const double x_max = 10.0;

    const double vx_min = -10.0 * sqrt(T_e / m_e);
    const double vx_max = 10.0 * sqrt(T_e / m_e);

    const double L = fabs(x_max - x_min);
    const double vx_width = fabs(vx_max - vx_min);

    //grids number
    const int grids_num = 200;
    const double grid_width = L / grids_num;

    //simulated steps & dt
    const int maxsteps = 100;
    const int time_ran = 0.0;
    const double timestep_condition = 0.1;

    //data path
    const string data_path = "./data/";

    //special parameters
    //for inhomogeneous initial space distribution
    const double uA_e = 0.6;
    const double uA_i = 0.4;
    const double u_k = 1.0;

    //electron/ion average number
    const double n_e_aver = N_e / L;
    const double n_i_aver = N_i / L;

    //calculated parameters
    const double q_e = -1.0 * sqrt(1.0 / n_e_aver);
    const double w_pe = sqrt(n_e_aver * q_e * q_e / m_e);
    const double w_pi = sqrt(n_i_aver * q_i * q_i / m_i);
    const double w_p = w_pe;
    const double dt = timestep_condition / w_p;
    const double lambda_D = sqrt(T_e / n_e_aver / q_e / q_e);

    void PrintParameters() const;

    double GetElecInitDistrib(double x,double v)
    {
        double ue = 1.0 + uA_e * cos( 2.0 * M_PI * x / L);
        double f = sqrt( m_e / (2 * M_PI * T_e) ) * exp(-m_e * (v * v) / 2 / T_e);
        return ue / L * f;
    }
    double max_probability_density = GetElecInitDistrib(0.0,0.0);

    double GetIonDensity(double x)
    {
        double ui = 1.0 + uA_i * cos( 2.0 * M_PI * x / L);
        return ui;
    }
};
#endif

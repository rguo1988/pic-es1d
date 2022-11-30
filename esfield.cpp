#include "esfield.h"
#include<iostream>

double LineInterpolation(double x, int k_gridpoint, double grid_width) //Linear interpolation
{
    double r = 0.0;
    if(fabs(x - k_gridpoint * grid_width) <= grid_width)
        r = (grid_width - fabs(x - k_gridpoint * grid_width)) / grid_width / grid_width;
    return r;
}

double QuadInterpolation(double x, int nearest_gridpoint, double grid_width, int index) //qs interpolation
{
    double  r = 0.0;
    switch(index)
    {
    case 0:
        r = pow(0.5 + (x - nearest_gridpoint * grid_width) / grid_width, 2) / grid_width / 2.0;
        break;
    case 1:
        r = (0.75 - pow((x - nearest_gridpoint * grid_width) / grid_width, 2)) / grid_width;
        break;
    case 2:
        r = pow(0.5 - (x - nearest_gridpoint * grid_width) / grid_width, 2) / grid_width / 2.0;
        break;
    }
    return r;
}

map<int, double> PartitionToGrid(double grid_width, int nx_grids, double x, unsigned int interpolation)
{
    map<int, double> partition;
    partition.clear();
    double pX = (x >= 0) ? x : x + grid_width * nx_grids;

    switch(interpolation)
    {
    case 1:
    {
        int low_x_grid = floor(pX / grid_width);
        int up_x_grid = (low_x_grid + 1 == nx_grids) ? 0 : (low_x_grid + 1);

        partition.insert(make_pair(low_x_grid, LineInterpolation(pX, low_x_grid, grid_width)));
        partition.insert(make_pair(up_x_grid,  LineInterpolation(pX, low_x_grid + 1, grid_width)));
    }
    break;

    case 2:
    {
        int k_grid = floor(pX / grid_width + 0.5);
        int nearest_gridpoint_x[3] = {k_grid - 1, k_grid, k_grid + 1};

        for(int i = 0; i < 3; i++)
        {
            int txNGP = nearest_gridpoint_x[i];
            //period condition
            if(txNGP >= nx_grids)
                txNGP -= nx_grids;
            if(txNGP < 0)
                txNGP += nx_grids;

            double contrib_to_grid_k = QuadInterpolation(pX, k_grid, grid_width, i);
            partition.insert(make_pair(txNGP, contrib_to_grid_k));
        }
    }
    break;

    }
    return partition;
}
// 3d array transform 1d

//int IndexFrom2dTo1d(int x, int y, int N) //N^3 is length of 1d array
//{
//return x + y * N;
//}

//int IndexFrom1dTo2d(int i, int N)
//{
//int x = i % N;
//int y = i / N;
//GridPoint temp(x, y);
//return temp;
//}
ElectricField::ElectricField(int _grids_num, double _grid_width):
    L(_grids_num * _grid_width), nx_grids(_grids_num), dx(_grid_width) {}

void ElectricField::SolveEFieldOnGrids(VectorXd charge)
{
    //charge density in grids and Fourier space
    fftw_complex *charge_dens_grid_fft;
    charge_dens_grid_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * pow(nx_grids, 1));
    fftw_complex *char_dens_in_fourier_space;
    char_dens_in_fourier_space = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * pow(nx_grids, 1));

    //electric potential in grids and Fourier space
    fftw_complex *elec_pot_grid_fft;//electric potential
    elec_pot_grid_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * pow(nx_grids, 1));
    fftw_complex *elec_pot_in_fourier_space;
    elec_pot_in_fourier_space = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * pow(nx_grids, 1));

    //input charge_dens_grid_fft
    for(int i = 0; i < pow(nx_grids, 1); i++)
    {
        charge_dens_grid_fft[i][0] = charge[i];
        charge_dens_grid_fft[i][1] = 0;
    }

    //perform fft
    fftw_plan p_forward;
    p_forward = fftw_plan_dft_1d(nx_grids, charge_dens_grid_fft, char_dens_in_fourier_space, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p_forward);//fft on charDensity
    //calculate electric potential

    elec_pot_in_fourier_space[0][0] = 0;
    elec_pot_in_fourier_space[0][1] = 0;

    for(int j = 1; j < nx_grids; j++)
    {
        double k = 2 * M_PI * j / L;
        double delta = k * dx / 2.0;
        double s = pow( sin(delta) / delta, 2);
        double as = k * k * s;
        elec_pot_in_fourier_space[j][0] = char_dens_in_fourier_space[j][0] /  as / pow(nx_grids, 1);
        elec_pot_in_fourier_space[j][1] = char_dens_in_fourier_space[j][1] /  as / pow(nx_grids, 1);
    }

    //perform inverse FFT to derive potential
    fftw_plan p_backward;
    p_backward = fftw_plan_dft_1d(nx_grids, elec_pot_in_fourier_space, elec_pot_grid_fft, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p_backward);

    phi.resize(nx_grids);
    Ex.resize(nx_grids);
    //transform 1d elec_pot_fft to elec_pot_grid
    for(int k = 0; k < pow(nx_grids, 1); k++)
    {
        //GridPoint gp = IndexFrom1dTo2d(k, nx_grids);
        //elec_pot_grid.insert(make_pair(gp, elec_pot_grid_fft[k][0]));
        phi[k] = elec_pot_grid_fft[k][0];
    }
    //calculate elecFieldGrids
    for(int i = 0; i < nx_grids; i++)
    {
        int i_minus = (i == 0) ? (nx_grids - 1) : (i - 1);
        int i_plus = (i == nx_grids - 1) ? 0 : (i + 1);
        Ex[i] = (phi[i_minus] - phi[i_plus]) / 2.0 / dx;
    }

    fftw_destroy_plan(p_forward);
    fftw_destroy_plan(p_backward);
    fftw_free(charge_dens_grid_fft);
    fftw_free(char_dens_in_fourier_space);
    fftw_free(elec_pot_grid_fft);
    fftw_free(elec_pot_in_fourier_space);
    fftw_cleanup();
}
double ElectricField::GetEx(int x_idx)
{
    return Ex[x_idx];
}
double ElectricField::GetPhi(int x_idx)
{
    return phi[x_idx];
}


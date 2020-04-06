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
//
//class GridPoint
GridPoint::GridPoint(int xx, int yy)
{
    x = xx;
    y = yy;
}

bool operator<(const GridPoint& lhs, const GridPoint& rhs)
{
    if(lhs.x < rhs.x)
        return true;
    else if(lhs.x >= rhs.x)
        return false;
    return false;
}//���ز������������ܵ�map�ؼ���;��С�߼���xyz->z*N*N+y*N+x*N

map<GridPoint, double> PartitionToGrid(double grid_width, int grids_num, PhaseSpace temp_ps, double xmin, unsigned int interpolation)
{
    map<GridPoint, double> partition;
    partition.clear();
    double pX = (temp_ps.x - xmin);
    if (pX < 0) pX += grid_width * grids_num;

    switch(interpolation)
    {
    case 1:
    {
        int low_x_grid = floor(pX / grid_width);
        int up_x_grid = low_x_grid + 1;

        vector<int> nearest_gridpoint_x = {low_x_grid, up_x_grid};
        for(auto xNGP : nearest_gridpoint_x)
        {
            int txNGP = xNGP;
            //period condition
            if(xNGP == grids_num)
                txNGP = 0;
            double contrib_to_grid_k = LineInterpolation(pX, xNGP, grid_width);

            //���䵽����
            GridPoint temp_gp(txNGP, 0);
            partition.insert(make_pair(temp_gp, contrib_to_grid_k));
        }

    }
    break;

    case 2:
    {
        int t_x_grid = floor(2 * pX / grid_width);
        int k_grid = (t_x_grid % 2 == 0) ? t_x_grid / 2 : (t_x_grid + 1) / 2;
        int nearest_gridpoint_x[3] = {k_grid - 1, k_grid, k_grid + 1};

        for(int i = 0; i < 3; i++)
        {
            int txNGP = nearest_gridpoint_x[i];
            //period condition
            if(txNGP >= grids_num)
                txNGP -= grids_num;
            if(txNGP < 0)
                txNGP += grids_num;

            //���䵽����
            double contrib_to_grid_k = QuadInterpolation(pX, k_grid, grid_width, i);
            GridPoint temp_gp(txNGP);
            partition.insert(make_pair(temp_gp, contrib_to_grid_k));
        }

    }
    break;

    }

    return partition;
}
// 3d array transform 1d

int IndexFrom2dTo1d(int x, int y, int N) //N^3 is length of 1d array
{
    return x + y * N;
}
GridPoint IndexFrom1dTo2d(int i, int N)
{
    int x = i % N;
    int y = i / N;
    GridPoint temp(x, y);
    return temp;
}
//class ElectricFieldGrid
/*
ElectricFieldGrids::ElectricFieldGrids()
{
    charge_dens_grid.clear();
    elec_pot_grid.clear();
    Ex_grid.clear();
}
*/

void ElectricFieldGrids::SetupEFieldOnGrids()
{
    //charge density in grids and Fourier space
    fftw_complex *charge_dens_grid_fft;
    charge_dens_grid_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * pow(grids_num, 1));
    fftw_complex *char_dens_in_fourier_space;
    char_dens_in_fourier_space = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * pow(grids_num, 1));

    //electric potential in grids and Fourier space
    fftw_complex *elec_pot_grid_fft;//electric potential
    elec_pot_grid_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * pow(grids_num, 1));
    fftw_complex *elec_pot_in_fourier_space;
    elec_pot_in_fourier_space = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * pow(grids_num, 1));

    //input charge_dens_grid_fft
    for(int i = 0; i < pow(grids_num, 1); i++)
    {
        charge_dens_grid_fft[i][0] = charge_dens_grid[IndexFrom1dTo2d(i, grids_num)];
        charge_dens_grid_fft[i][1] = 0;
    }

    //perform fft
    fftw_plan p_forward;
    p_forward = fftw_plan_dft_1d(grids_num, charge_dens_grid_fft, char_dens_in_fourier_space, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p_forward);//fft on charDensity
    //calculate electric potential

    elec_pot_in_fourier_space[0][0] = 0;
    elec_pot_in_fourier_space[0][1] = 0;

    for(int j = 1; j < grids_num; j++)
    {
        //double kx = 2 * M_PI * j / (double)grids_num;
        double k = 2 * M_PI * j / L;
        double delta = k * grid_width / 2.0;
        double s = pow( sin(delta) / delta, 2);
        double as = k * k * s;
        elec_pot_in_fourier_space[j][0] = char_dens_in_fourier_space[j][0] /  as / pow(grids_num, 1);
        elec_pot_in_fourier_space[j][1] = char_dens_in_fourier_space[j][1] /  as / pow(grids_num, 1);
    }

    //perform inverse FFT to derive potential
    fftw_plan p_backward;
    p_backward = fftw_plan_dft_1d(grids_num, elec_pot_in_fourier_space, elec_pot_grid_fft, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p_backward);

    elec_pot_grid.clear();
    Ex_grid.clear();
    //transform 1d elec_pot_fft to elec_pot_grid
    for(int k = 0; k < pow(grids_num, 1); k++)
    {
        GridPoint gp = IndexFrom1dTo2d(k, grids_num);
        elec_pot_grid.insert(make_pair(gp, elec_pot_grid_fft[k][0]));
    }
    //calculate elecFieldGrids
    for(auto p : elec_pot_grid)
    {
        int x = p.first.x;

        int x_1 = (p.first.x == 0) ? (grids_num - 1) : (x - 1);

        int x1 = (p.first.x == grids_num - 1) ? 0 : (x + 1);

        GridPoint xy(x);
        GridPoint x1y(x1);
        GridPoint x_1y(x_1);

        double tempEx;
        tempEx = (elec_pot_grid[x_1y] - elec_pot_grid[x1y]) / 2 / grid_width;

        Ex_grid.insert(make_pair(xy, tempEx));
    }

    fftw_destroy_plan(p_forward);
    fftw_destroy_plan(p_backward);
    fftw_free(charge_dens_grid_fft);
    fftw_free(char_dens_in_fourier_space);
    fftw_free(elec_pot_grid_fft);
    fftw_free(elec_pot_in_fourier_space);
    fftw_cleanup();
}
double ElectricFieldGrids::AppliedElecFields(double x, double t, double amplitude, double k, double w)
{
    double r = amplitude * cos(k * x - w * t);
    return r;
}

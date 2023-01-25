#include "poisson_solver.h"

PoissonSolverPeriodicBC_FFTW::PoissonSolverPeriodicBC_FFTW(int _nx, double _grid_width):
    nx_grids(_nx - 1), dx(_grid_width), L((_nx - 1) * _grid_width)
{

    phi.resize(nx_grids);
    Ex.resize(nx_grids);
}

void PoissonSolverPeriodicBC_FFTW::Solve(VectorXd charge)
{
    //charge density in grids and Fourier space
    fftw_complex *charge_dens_grid_fft;
    charge_dens_grid_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx_grids);
    fftw_complex *char_dens_in_fourier_space;
    char_dens_in_fourier_space = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx_grids);

    //electric potential in grids and Fourier space
    fftw_complex *elec_pot_grid_fft;//electric potential
    elec_pot_grid_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx_grids);
    fftw_complex *elec_pot_in_fourier_space;
    elec_pot_in_fourier_space = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx_grids);

    //input charge_dens_grid_fft
    for(int i = 0; i < nx_grids; i++)
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
        elec_pot_in_fourier_space[j][0] = char_dens_in_fourier_space[j][0] /  as / nx_grids;
        elec_pot_in_fourier_space[j][1] = char_dens_in_fourier_space[j][1] /  as / nx_grids;
    }

    //perform inverse FFT to derive potential
    fftw_plan p_backward;
    p_backward = fftw_plan_dft_1d(nx_grids, elec_pot_in_fourier_space, elec_pot_grid_fft, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p_backward);

    //transform 1d elec_pot_fft to elec_pot_grid
    for(int k = 0; k < nx_grids; k++)
    {
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
double PoissonSolverPeriodicBC_FFTW::GetE(int x_idx)
{
    return Ex[x_idx];
}
double PoissonSolverPeriodicBC_FFTW::GetPhi(int x_idx)
{
    return phi[x_idx];
}

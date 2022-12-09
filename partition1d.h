#ifndef _partition1d_h
#define _partition1d_h

#include <cmath>
#include <map>

using namespace std;

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
#endif

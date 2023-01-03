//general head file
#include <iostream>
#include<omp.h>

//project head file
#include"plasma.h"

using namespace std;

int main()
{
    double start_time = omp_get_wtime();
    //creat plasma
    PlasmaSystem plasma;

    plasma.Run();
    double stop_time = omp_get_wtime();
    cout << "  Using Time: " << stop_time - start_time << endl;
    return 0;
}

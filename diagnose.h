//测试，调试，以及输出文件
#ifndef _diagnose_h
#define _diagnose_h
#include <fstream>
#include <string>
#include "esfield.h"

using namespace std;

double OutputTotCharge(ElectricField cdg);

void OutputData(string filename, vector<double> a);
vector<double> GetParticlesVX(Particles testp);
vector<double> GetParticlesX(Particles testp);

void OutputNetCharge(ElectricField efg);

void OutputChargeGrids(string filename, ElectricField cdg);
void OutputElectricPotentialGrids(string filename, ElectricField efg);
void OutputElectricFieldGrids(string filename, ElectricField efg);

#endif

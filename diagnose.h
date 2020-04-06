//测试，调试，以及输出文件
#ifndef _diagnose_h
#define _diagnose_h
#include <fstream>
#include <string>
#include "esfield.h"

using namespace std;

double OutputTotCharge(ElectricFieldGrids cdg);

void OutputData(string filename, vector<double> a);
vector<double> GetParticlesVX(Particles testp);
vector<double> GetParticlesX(Particles testp);

void OutputNetCharge(ElectricFieldGrids efg);

void OutputChargeGrids(string filename, ElectricFieldGrids cdg);
void OutputElectricPotentialGrids(string filename, ElectricFieldGrids efg);
void OutputElectricFieldGrids(string filename, ElectricFieldGrids efg);

#endif

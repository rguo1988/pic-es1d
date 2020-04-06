//������Ķ���
//
#include "particles.h"

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

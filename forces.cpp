#include "forces.hpp"

#include <cmath>
#include <raylib.h>

#include "system.hpp"

Gravity::Gravity(double a_)
    : a(a_)
{
}

void Gravity::Apply(System& system)
{
    for (int i = 0; i < system.N; i++)
        system.force.At(i*system.DF+1) += (1.0/system.massInv.At(i)) * a;
}

Spring::Spring(int a_, int b_, double len_, double k_)
    : a(a_), b(b_), len(len_), k(k_)
{
}

void Spring::Apply(System& system)
{
    double ax = system.pos.At(a*system.DF);
    double ay = system.pos.At(a*system.DF+1);

    double bx = system.pos.At(b*system.DF);
    double by = system.pos.At(b*system.DF+1);

    double dx = (ax - bx);
    double dy = (ay - by);

    double d = std::sqrt(dx*dx + dy*dy);

    double x = d - len;

    double F = -k * x;

    system.force.At(a*system.DF) += -1.0*dx/d*F;
    system.force.At(a*system.DF+1) += -1.0*dy/d*F;

    system.force.At(b*system.DF) += dx/d*F;
    system.force.At(b*system.DF+1) += dy/d*F;
}

void Spring::Draw(System& system)
{
    double ax = system.pos.At(a*system.DF);
    double ay = system.pos.At(a*system.DF+1);

    double bx = system.pos.At(b*system.DF);
    double by = system.pos.At(b*system.DF+1);

    DrawLine(ax, ay, bx, by, YELLOW);
}

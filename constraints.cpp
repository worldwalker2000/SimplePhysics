#include "constraints.hpp"

#include <cmath>

#include "system.hpp"

PositionConstraint::PositionConstraint(int a_, double x_, double y_)
    : a(a_), x(x_), y(y_)
{
}

// constraint value
void PositionConstraint::C(System& system, int i)
{
    double ax = system.pos.At(a*system.DF);
    double ay = system.pos.At(a*system.DF+1);

    double dx = x - ax;
    double dy = y - ay;

    double d = std::sqrt(dx*dx + dy*dy);

    system.C.At(i) = d;
}

// time derivative of C
void PositionConstraint::Cd(System& system, int i)
{
    double ax = system.pos.At(a*system.DF);
    double ay = system.pos.At(a*system.DF+1);

    double vx = system.vel.At(a*system.DF);
    double vy = system.vel.At(a*system.DF+1);

    double dx = x - ax;
    double dy = y - ay;

    double d = std::sqrt(dx*dx + dy*dy);

    double p1 = (x - ax) * vx;
    double p2 = (y - ay) * vy;

    if (d == 0.0) system.Cd.At(i) = 0.0;
    else system.Cd.At(i) = (-1.0*(p1+p2))/d;
}

// derivative of C wrt pos (q)
void PositionConstraint::J(System& system, int i)
{
    double ax = system.pos.At(a*system.DF);
    double ay = system.pos.At(a*system.DF+1);

    double dx = x - ax;
    double dy = y - ay;

    double d = std::sqrt(dx*dx + dy*dy);

    if (d == 0.0)
    {
        system.J.At(i, a*system.DF) = 0.0;
        system.J.At(i, a*system.DF+1) = 0.0;
    }
    else
    {
        system.J.At(i, a*system.DF) = (-1.0*(x-ax))/d;
        system.J.At(i, a*system.DF+1) = (-1.0*(y-ay))/d;
    }
}

// derivative of Cd wrt pos (q) or time derivative of J
void PositionConstraint::Jd(System& system, int i)
{
    double ax = system.pos.At(a*system.DF);
    double ay = system.pos.At(a*system.DF+1);

    double vx = system.vel.At(a*system.DF);
    double vy = system.vel.At(a*system.DF+1);

    double dx = x - ax;
    double dy = y - ay;

    double d = std::sqrt(dx*dx + dy*dy);

    double p1 = (x - ax) * vx;
    double p2 = (y - ay) * vy;

    double Csq = system.C.At(i)*system.C.At(i);

    if (Csq == 0.0)
    {
        system.Jd.At(i, a*system.DF) = 0.0;
        system.Jd.At(i, a*system.DF+1) = 0.0;
    }
    else
    {
        system.Jd.At(i, a*system.DF) = ((system.C.At(i)*vx)-((-1.0*(p1+p2))*system.J.At(i, a*system.DF)))/Csq;
        system.Jd.At(i, a*system.DF+1) = ((system.C.At(i)*vy)-((-1.0*(p1+p2))*system.J.At(i, a*system.DF+1)))/Csq;
    }
}


DistanceConstraint::DistanceConstraint(int a_, int b_, double dist_)
    : a(a_), b(b_), dist(dist_)
{
}

// constraint value
void DistanceConstraint::C(System& system, int i)
{
    double ax = system.pos.At(a*system.DF);
    double ay = system.pos.At(a*system.DF+1);

    double bx = system.pos.At(b*system.DF);
    double by = system.pos.At(b*system.DF+1);

    double dx = bx - ax;
    double dy = by - ay;

    double d = std::sqrt(dx*dx + dy*dy);

    // TODO: which way should this subtraction be?
    system.C.At(i) = d - dist;
}

// time derivative of C
void DistanceConstraint::Cd(System& system, int i)
{
    double ax = system.pos.At(a*system.DF);
    double ay = system.pos.At(a*system.DF+1);

    double avx = system.vel.At(a*system.DF);
    double avy = system.vel.At(a*system.DF+1);

    double bx = system.pos.At(b*system.DF);
    double by = system.pos.At(b*system.DF+1);

    double bvx = system.vel.At(b*system.DF);
    double bvy = system.vel.At(b*system.DF+1);

    double dx = bx - ax;
    double dy = by - ay;

    double d = std::sqrt(dx*dx + dy*dy);

    double top = (bx-ax)*(bvx-avx)+(by-ay)*(bvy-avy);

    // TODO: which way should this subtraction be?
    if (d == 0.0)
    {
        system.Cd.At(i) = 0.0;
    }
    else
    {
        system.Cd.At(i) = ((top)/d);
    }
}

// derivative of C wrt system.pos (q)
void DistanceConstraint::J(System& system, int i)
{
    double ax = system.pos.At(a*system.DF);
    double ay = system.pos.At(a*system.DF+1);

    double bx = system.pos.At(b*system.DF);
    double by = system.pos.At(b*system.DF+1);

    double dx = bx - ax;
    double dy = by - ay;

    double d = std::sqrt(dx*dx + dy*dy);

    if (d == 0.0)
    {
        system.J.At(i, a*system.DF) = 0.0;
        system.J.At(i, a*system.DF+1) = 0.0;

        system.J.At(i, b*system.DF) = 0.0;
        system.J.At(i, b*system.DF+1) = 0.0;
    }
    else
    {
        system.J.At(i, a*system.DF) = -1.0*dx/d;
        system.J.At(i, a*system.DF+1) = -1.0*dy/d;

        system.J.At(i, b*system.DF) = dx/d;
        system.J.At(i, b*system.DF+1) = dy/d;
    }
}

// derivative of Cd wrt system.pos (q) or time derivative of J
void DistanceConstraint::Jd(System& system, int i)
{
    double ax = system.pos.At(a*system.DF);
    double ay = system.pos.At(a*system.DF+1);

    double avx = system.vel.At(a*system.DF);
    double avy = system.vel.At(a*system.DF+1);

    double bx = system.pos.At(b*system.DF);
    double by = system.pos.At(b*system.DF+1);

    double bvx = system.vel.At(b*system.DF);
    double bvy = system.vel.At(b*system.DF+1);

    double dx = bx - ax;
    double dy = by - ay;

    double d = std::sqrt(dx*dx + dy*dy);
    double top = (bx-ax)*(bvx-avx)+(by-ay)*(bvy-avy);

    double dsq = d*d;

    if (dsq == 0.0)
    {
        system.Jd.At(i, a*system.DF) = 0.0;
        system.Jd.At(i, a*system.DF+1) = 0.0;

        system.Jd.At(i, b*system.DF) = 0.0;
        system.Jd.At(i, b*system.DF+1) = 0.0;
    }
    else
    {
        system.Jd.At(i, a*system.DF) = (-1.0*d*(bvx-avx)-top*system.J.At(i, a*system.DF))/dsq;
        system.Jd.At(i, a*system.DF+1) = (-1.0*d*(bvy-avy)-top*system.J.At(i, a*system.DF+1))/dsq;

        system.Jd.At(i, b*system.DF) = (d*(bvx-avx)-top*system.J.At(i, b*system.DF))/dsq;
        system.Jd.At(i, b*system.DF+1) = (d*(bvy-avy)-top*system.J.At(i, b*system.DF+1))/dsq;
    }
}

void DistanceConstraint::Draw(System& system)
{
    double ax = system.pos.At(a*system.DF);
    double ay = system.pos.At(a*system.DF+1);

    double bx = system.pos.At(b*system.DF);
    double by = system.pos.At(b*system.DF+1);

    DrawLine(ax, ay, bx, by, BLUE);
}

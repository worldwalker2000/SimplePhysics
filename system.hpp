#pragma once

#include <vector>

#include "la.hpp"

#include "forces.hpp"
#include "constraints.hpp"

struct Particle
{
    double x;
    double y;
    double m;
};

struct System
{
    const int N;
    const int DF;
    const int NC;
    const int NF;

    Vec pos;
    Vec vel;
    Vec force;
    Vec massInv;
    Mat W;

    Vec C;
    Vec Cd;
    Mat J;
    Mat Jd;

    const double ks;
    const double kd;

    std::vector<Force*> forces;
    std::vector<Constraint*> constraints;

    double totalError;

    System(const std::vector<Particle>& particles, std::vector<Force*> forces_, std::vector<Constraint*> constraints_);
    ~System();

    void ResetConstraints();

    void Step(double dt, int steps);

    void Draw();
};

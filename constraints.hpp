#pragma once

#include <raylib.h>

#include "la.hpp"

class System;

struct Constraint
{
    virtual ~Constraint() = default;

    virtual void C(System& system, int i) = 0;
    virtual void Cd(System& system, int i) = 0;

    virtual void J(System& system, int i) = 0;
    virtual void Jd(System& system, int i) = 0;

    virtual void Draw(System& system) {  };
};

struct PositionConstraint : public Constraint
{
    int a;
    double x, y;

    PositionConstraint(int a_, double x_, double y_);

    virtual void C(System& system, int i) override;
    virtual void Cd(System& system, int i) override;
    virtual void J(System& system, int i) override;
    virtual void Jd(System& system, int i) override;
};

struct DistanceConstraint : public Constraint
{
    int a, b;
    double dist;

    DistanceConstraint(int a_, int b_, double dist_);

    virtual void C(System& system, int i) override;
    virtual void Cd(System& system, int i) override;
    virtual void J(System& system, int i) override;
    virtual void Jd(System& system, int i) override;

    virtual void Draw(System& system) override;
};

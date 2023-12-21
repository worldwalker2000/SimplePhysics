#pragma once

#include "la.hpp"

class System;

struct Force
{
    virtual ~Force() = default;

    virtual void Apply(System& system) = 0;
    virtual void Draw(System& system) {  }
};

struct Gravity : public Force
{
    double a;

    Gravity(double a_);

    virtual void Apply(System& system) override;
};

struct Spring : public Force
{
    int a, b;

    double len;
    double k;

    Spring(int a_, int b_, double len_, double k_);

    virtual void Apply(System& system) override;
    virtual void Draw(System& system) override;
};

#include "system.hpp"

#include <cmath>
#include <raylib.h>

System::System(const std::vector<Particle>& particles, std::vector<Force*> forces_, std::vector<Constraint*> constraints_)
    : N(particles.size())
    , DF(2)
    , NC(constraints_.size())
    , NF(forces_.size())
    , pos(N*DF)
    , vel(N*DF)
    , force(N*DF)
    , massInv(N*DF)
    , W(N*DF, N*DF)
    , C(NC)
    , Cd(NC)
    , J(NC, N*DF)
    , Jd(NC, N*DF)
    , ks(0.1)
    , kd(0.1)
    , forces(forces_)
    , constraints(constraints_)
    , totalError(0.0)
{
    for (int i = 0; i < N; i++)
    {
        pos.At(i*DF) = particles[i].x;
        pos.At(i*DF+1) = particles[i].y;

        vel.At(i*DF) = 0.0;
        vel.At(i*DF+1) = 0.0;

        massInv.At(i*DF) = massInv.At(i*DF+1) = 1.0/particles[i].m;

        W.At(i*DF, i*DF) = W.At(i*DF+1, i*DF+1) = 1.0/particles[i].m;
    }

    ResetConstraints();
}

System::~System()
{
    /*
    for (int i = 0; i < NF; i++)
        delete forces[i];

    for (int i = 0; i < NC; i++)
        delete constraints[i];
        */
}

void System::ResetConstraints()
{
    C.Zero(); Cd.Zero(); J.Zero(); Jd.Zero();
}

void System::Step(double dt, int steps)
{
    for (int step = 0; step < steps; step++)
    {
        force.Zero();

        for (int i = 0; i < NF; i++)
            forces[i]->Apply(*this);

        // Constraints
        {
            ResetConstraints();

            for (int i = 0; i < NC; i++)
            {
                constraints[i]->C(*this, i);
                constraints[i]->Cd(*this, i);
                constraints[i]->J(*this, i);
                constraints[i]->Jd(*this, i);
            }

            // (J*W*Jt) * l = -Jd*qd - J*W*Q - ks*C - kd * Cd

            // J*W*Jt
            Mat Jt(J); Jt.Transpose();
            Mat A = J*W*Jt;

            // Jd*qd
            Vec Jdqd = Jd*vel;

            // J*W*Q
            Vec JWQ = J*W*force;

            // -Jdqd - JWQ - ks*C - kd*Cd
            Vec b = -1.0*Jdqd + -1.0*JWQ + -1.0*(ks*C) + -1.0*(kd*Cd);

            // Solve A*l=b
            Vec l = Mat::Solve(A, b);

            // Qh = Jt*l
            Vec Qh = Jt*l;

            // force + Qh
            force += Qh;
        }

        Vec acl = force*massInv;

        // TODO: rk4 integrator

        /*
            float rk4(float x0, float y0, float x1, float step)
            {
                    // k1 = hf(x0, y0)
                    // k2 = hf[x0 + (0.5)h, y0 + (0.5)k1]
                    // k3 = hf[x0 + (0.5)h, y0 + (0.5)k2]
                    // k4 = hf(x0 + h, y0 + k3)

                    // y1 = y0 + (1/6) (k1 + 2k2 + 2k3 + k4)

                float y = y0;

                for (float a = x0; a < x1; a += step) {
                    float k1 = step * dydx(a, y);
                    float k2 = step * dydx(a + 0.5f*step, y + 0.5f*k1);
                    float k3 = step * dydx(a + 0.5f*step, y + 0.5f*k2);
                    float k4 = step * dydx(a + step, y + k3);

                    y = y + (1.0f/6.0f) * (k1+2*k2+2*k3+k4);
                }


                return y;
            }
        */

        vel += acl*(dt/steps);
        pos += vel*(dt/steps);

        totalError = 0.0;
        for (int i = 0; i < NC; i++)
            totalError += std::abs(C.At(i));
    }
}

void System::Draw()
{
    for (int i = 0; i < N; i++) DrawCircle(pos.At(i*DF), pos.At(i*DF+1), 20, RED);

    for (int i = 0; i < NF; i++)
        forces[i]->Draw(*this);

    for (int i = 0; i < NC; i++)
        constraints[i]->Draw(*this);
}

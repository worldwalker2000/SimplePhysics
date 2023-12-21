#include <cmath>
#include <iostream>
#include <vector>
#include <chrono>
#include <thread>
#include <cstdlib>
#include <cstring>
#include <glm/glm.hpp>
#include <raylib.h>

#include "constraints.hpp"
#include "forces.hpp"
#include "la.hpp"
#include "system.hpp"

typedef std::chrono::high_resolution_clock Clock;

constexpr int width = 1300;
constexpr int height = 800;

double Dist(double x, double y, double cx, double cy)
{
    double dx = x - cx;
    double dy = y - cy;

    return std::sqrt(dx*dx + dy*dy);
}

bool InCircle(double x, double y, double cx, double cy, double r)
{
    return Dist(x, y, cx, cy) < r;
}

enum State
{
    BUILD = 0,
    SIM,
};

enum Tool
{
    MASS1 = 1,
    MASS3 = 3,
    MASS5 = 5,
    MASS10 = 10,
    SPRING15 = 15,
    SPRING100 = 100,
    SPRING500 = 500,
    SPRING1000 = 1000,
    ROD = 2000,
    HOLD = 3000,
};

int main()
{
    SetTraceLogLevel(LOG_NONE);
    InitWindow(width, height, "Simple Physics");
    // SetWindowMonitor(0);
    SetTargetFPS(60);

    Clock::time_point lastTime = Clock::now();
    double timeAcc = 0.0f;

    constexpr int dist = 200;

    std::vector<Particle> particles;

    std::vector<Force*> forces = {
        new Gravity(200.0f),
    };

    std::vector<Constraint*> constraints;

    System* system = NULL;

    State state = BUILD;
    Tool tool = MASS1;
    int a = -1;
    int b = -1;

    while (!WindowShouldClose())
    {
        std::chrono::high_resolution_clock::time_point current = std::chrono::high_resolution_clock::now();
        long long deltaUs = std::chrono::duration_cast<std::chrono::microseconds>(current - lastTime).count();
        double deltaMs = deltaUs / 1000.0f;
        double dt = std::chrono::duration_cast<std::chrono::microseconds>(current - lastTime).count() / 1000000.0f;
        double fps = 1 / dt;
        lastTime = current;
        timeAcc += dt;

        constexpr int steps = 1000;

        if (timeAcc >= 1.0) { std::printf("fps: %0.2f\n", fps); timeAcc = 0.0; }

        switch (state)
        {
            case State::BUILD:
                {
                    if (system) { delete system; system = NULL; }

                    double x = GetMouseX();
                    double y = GetMouseY();
                    if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT))
                    {
                        switch (tool)
                        {
                            case MASS1:
                            case MASS3:
                            case MASS5:
                            case MASS10:
                                particles.push_back({ .x = x, .y = y, .m = (double) tool });
                                break;
                            case SPRING15:
                            case SPRING100:
                            case SPRING500:
                            case SPRING1000:
                                if (a == -1)
                                {
                                    for (int i = 0; i < particles.size(); i++)
                                        if (InCircle(x, y, particles[i].x, particles[i].y, 20)) a = i;
                                }
                                else
                                {
                                    for (int i = 0; i < particles.size(); i++)
                                        if (InCircle(x, y, particles[i].x, particles[i].y, 20)) b = i;

                                    if (b != -1)
                                    {
                                        forces.push_back(new Spring(a, b, Dist(particles[a].x, particles[a].y, particles[b].x, particles[b].y), -1.0 * (double) tool));
                                        a = b = -1;
                                    }
                                }
                                break;
                            case ROD:
                                if (a == -1)
                                {
                                    for (int i = 0; i < particles.size(); i++)
                                        if (InCircle(x, y, particles[i].x, particles[i].y, 20)) a = i;
                                }
                                else
                                {
                                    for (int i = 0; i < particles.size(); i++)
                                        if (InCircle(x, y, particles[i].x, particles[i].y, 20)) b = i;

                                    if (b != -1)
                                    {
                                        constraints.push_back(new DistanceConstraint(a, b, Dist(particles[a].x, particles[a].y, particles[b].x, particles[b].y)));
                                        a = b = -1;
                                    }
                                }
                                break;
                            case HOLD:
                                {
                                    int p = -1;
                                    for (int i = 0; i < particles.size(); i++)
                                        if (InCircle(x, y, particles[i].x, particles[i].y, 20)) p = i;

                                    if (p != -1)
                                        constraints.push_back(new PositionConstraint(p, particles[p].x, particles[p].y));
                                }
                                break;
                        }
                    }

                    if (IsKeyPressed(KEY_ONE)) tool = MASS1;
                    if (IsKeyPressed(KEY_TWO)) tool = MASS3;
                    if (IsKeyPressed(KEY_THREE)) tool = MASS5;
                    if (IsKeyPressed(KEY_FOUR)) tool = MASS10;
                    if (IsKeyPressed(KEY_FIVE)) tool = SPRING15;
                    if (IsKeyPressed(KEY_SIX)) tool = SPRING100;
                    if (IsKeyPressed(KEY_SEVEN)) tool = SPRING500;
                    if (IsKeyPressed(KEY_EIGHT)) tool = SPRING1000;
                    if (IsKeyPressed(KEY_NINE)) tool = ROD;
                    if (IsKeyPressed(KEY_ZERO)) tool = HOLD;

                    if (IsKeyPressed(KEY_C)) a = b = -1;

                    if (IsKeyPressed(KEY_SPACE)) state = SIM;

                    BeginDrawing();
                    {
                        ClearBackground(BLACK);

                        if (a != -1) DrawText("(*)", width - 30, 10, 20, WHITE);

                        const char* text = "Use keys to selet tools:\n"
                            "1 = MASS 1\n"
                            "2 = MASS 3\n"
                            "3 = MASS 5\n"
                            "4 = MASS 10\n"
                            "5 = SPRING 15\n"
                            "6 = SPRING 100\n"
                            "7 = SPRING 500\n"
                            "8 = SPRING 1000\n"
                            "9 = ROD\n"
                            "0 = HOLD\n";
                        DrawText(text, 10, 10, 20, WHITE);

                        for (int i = 0; i < particles.size(); i++) DrawCircle(particles[i].x, particles[i].y, 20, RED);

                        for (int i = 0; i < forces.size(); i++)
                        {
                            Spring* s;
                            if ((s = dynamic_cast<Spring*>(forces[i])))
                                DrawLine(particles[s->a].x, particles[s->a].y, particles[s->b].x, particles[s->b].y, YELLOW);
                        }

                        for (int i = 0; i < constraints.size(); i++)
                        {
                            PositionConstraint* p;
                            if ((p = dynamic_cast<PositionConstraint*>(constraints[i])))
                                DrawCircle(particles[p->a].x, particles[p->a].y, 20, ORANGE);

                            DistanceConstraint* d;
                            if ((d = dynamic_cast<DistanceConstraint*>(constraints[i])))
                                DrawLine(particles[d->a].x, particles[d->a].y, particles[d->b].x, particles[d->b].y, BLUE);
                        }

                    }
                    EndDrawing();
                }
                break;
            case State::SIM:
                {
                    if (!system) system = new System(particles, forces, constraints);

                    if (IsKeyPressed(KEY_SPACE)) state = BUILD;

                    system->Step(dt, 10000);

                    BeginDrawing();
                    {
                        ClearBackground(BLACK);

                        system->Draw();
                    }
                    EndDrawing();
                }
                break;
        }

        // std::this_thread::sleep_for(std::chrono::milliseconds(3));
    }

    if (system) delete system;

    for (int i = 0; i < forces.size(); i++) delete forces[i];
    for (int i = 0; i < constraints.size(); i++) delete constraints[i];

    CloseWindow();

    return 0;
}

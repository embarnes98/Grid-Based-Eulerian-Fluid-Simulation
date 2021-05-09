#ifndef FLUID_SIM_SFML_H_
#define FLUID_SIM_SFML_H_
#define Nx 256
#define Ny 256
#define ITER 4
#define DENSITY_TOLERANCE 0.5f
#define SPEED_TOLERANCE 0.7f
#define IN_CIRCLE(x, y, xc, yc, r) (x - xc) * (x - xc) + (y - yc) * (y - yc) < r * r
#define IX(_x, _y) ((_x) + ((_y)*(Nx)))
#define IFATSURFACE if (!(m_optimised) || m_surface[IX(i, j)]) {
#define FOR_EACH_CELL for (j = m_topEdge; j <= m_bottomEdge; j++) {for (i = m_leftEdge; i <= m_rightEdge; i++) {IFATSURFACE
#define END_FOR }}}
#include <cstdint>
#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <assert.h> 

class FluidGrid
{
  public:
    FluidGrid(float _diffusion, float _viscosity);
    ~FluidGrid();
    void step(float _dt);
    void addDensity(int _x, int _y, float _densityDelta);
    void addVelocity(int _x, int _y, float _deltavx, float _deltavy);
    int getDensity(int _i);
    void updateLocalGrid(int _x, int _y);
    int getSurfaceFlag(int _x, int _y);
    int getSurfaceFlag(int _i);
    int getLeftEdge();
    int getRightEdge();
    int getTopEdge();
    int getBottomEdge();
    int isOptimised();
  private:
    float m_diffusion;
    float m_viscosity;
    float m_s[Nx * Ny];
    float m_density[Nx * Ny];
    float m_vx[Nx * Ny];
    float m_vy[Nx * Ny];
    float m_vx0[Nx * Ny];
    float m_vy0[Nx * Ny];
    int m_leftEdge;
    int m_rightEdge;
    int m_topEdge;
    int m_bottomEdge;
    int m_surface[Nx * Ny];
    int m_optimised = 1;
    void prune();
    void setBoundary(uint8_t _b, float* _x);
    void linearSolve(uint8_t _b, float* _x, float* _x0, float _a, float _c);
    void diffuse(uint8_t _b, float* _x, float* _x0, float _diff, float _dt);
    void project(float* _velocX, float* _velocY, float* _p, float* _div);
    void advect(uint8_t _b, float* _d, float* _d0,  float* _velocX, float* _velocY, float _dt);
};

#endif

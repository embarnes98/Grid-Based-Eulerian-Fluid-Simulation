#include "Fluid-Sim-SFML.h"


FluidGrid::FluidGrid(float _diffusion, float _viscosity)
{
  std::cout << "Creating FluidGrid...\n";
  m_diffusion = _diffusion;
  m_viscosity = _viscosity;
  auto numPixels = Nx * Ny;
  for (auto i = 0; i < numPixels; i++)
  {
    m_s[i] = 0.f;
    m_density[i] = 0.f;
    m_vx[i] = 0.f;
    m_vy[i] = 0.f;
    m_vx0[i] = 0.f;
    m_vy0[i] = 0.f;
    m_surface[i] = 0;
  }
  if (m_optimised)
  {
    m_leftEdge = Nx/2;
    m_rightEdge = Nx/2;
    m_topEdge = Ny/2;
    m_bottomEdge = Ny/2;
  }
  else
  {
    m_leftEdge = 1;
    m_rightEdge = Nx - 2;
    m_topEdge = 1;
    m_bottomEdge = Ny - 2;
  }
}



FluidGrid::~FluidGrid()
{
  std::cout << "Destroying FluidGrid.\n";
}

void FluidGrid::addDensity(int _x, int _y, float _densityDelta)
{
  m_density[IX(_x, _y)] += _densityDelta;
}

void FluidGrid::addVelocity(int _x, int _y, float _deltavx, float _deltavy)
{
  auto i = IX(_x, _y);
  m_vx[i] += _deltavx;
  m_vy[i] += _deltavy;
}

// Seems to be working as expected
void FluidGrid::updateLocalGrid(int _x, int _y)
{
  auto spr = 4;
  m_leftEdge = std::max(std::min(m_leftEdge, _x - spr), 1);
  m_rightEdge = std::min(std::max(m_rightEdge, _x + spr), Nx - 2);
  m_topEdge = std::max(std::min(m_topEdge, _y - spr), 1);
  m_bottomEdge = std::min(std::max(m_bottomEdge, _y + spr), Ny - 2);
}

void FluidGrid::prune()
{ 
  int i, j;
  bool edgeExtended;
  auto tolerance = 1.f;
  auto halfTolerance = tolerance/2.f;
  float abs_vx, abs_vy;

  // Simulation edges calculation

  if (m_leftEdge > 0)
  {
    edgeExtended = false;
    for (i = m_topEdge; i <= m_bottomEdge; i++)
    {
      if (m_density[IX(m_leftEdge + 1, i)] > tolerance || m_density[IX(m_leftEdge + 2, i)] > tolerance)
      {
        m_leftEdge = std::max(m_leftEdge - 1, 1);
        edgeExtended = true;
        break;
      }
    }
    if (!edgeExtended)
    {
      m_leftEdge = std::min(Nx/2, m_leftEdge + 1);
    }
  }
  if (m_rightEdge < Nx - 1)
  {
    edgeExtended = false;
    for (i = m_topEdge; i <= m_bottomEdge; i++)
    {
      if (m_density[IX(m_rightEdge - 1, i)] > tolerance || m_density[IX(m_rightEdge - 2, i)] > tolerance) 
      {
        m_rightEdge = std::min(m_rightEdge + 1, Nx - 2); // changed from N to Nx - 1
        edgeExtended = true;
        break;
      }
    }
    if (!edgeExtended)
    {
      m_rightEdge = std::max(Nx/2, m_rightEdge - 1);
    }
  }
  if (m_topEdge > 0)
  {
    edgeExtended = false;
    for (i = m_leftEdge; i <= m_rightEdge; i++)
    {
      if (m_density[IX(i, m_topEdge + 1)] > tolerance || m_density[IX(i, m_topEdge + 2)] > tolerance) 
      {
        m_topEdge = std::max(m_topEdge - 1, 1);
        edgeExtended = true;
        break;
      }
    }
    if (!edgeExtended)
    {
      m_topEdge = std::min(Ny/2, m_topEdge + 1);
    }
  }
  if (m_bottomEdge < Ny - 1)
  {
    edgeExtended = false;
    for (i = m_leftEdge; i <= m_rightEdge; i++)
    {
      if (m_density[IX(i, m_bottomEdge - 1)] > 0 || m_density[IX(i, m_bottomEdge - 2)] > 0 )
      {
        m_bottomEdge = std::min(m_bottomEdge + 1, Ny - 2);
        edgeExtended = true;
        break;
      }
    }
    if (!edgeExtended)
    {
      m_bottomEdge = std::max(Ny/2, m_bottomEdge - 1);
    }
  }

  // Simulation surface calculation

  for (j = m_topEdge; j <= m_bottomEdge; j++) 
  {
    for (i = m_leftEdge; i <= m_rightEdge; i++) 
    {
      // if the density of this cell and its neighbours is small
      // enough, set its density and velocity to zero and don't
      // simulate it
      if  (
            m_density[IX(i, j    )] < tolerance &&
            m_density[IX(i, j + 1)] < tolerance &&
            m_density[IX(i, j - 1)] < tolerance &&
            m_density[IX(i + 1, j)] < tolerance &&
            m_density[IX(i - 1, j)] < tolerance 
          )
      { 
        m_surface[IX(i, j)] = 0;
        m_density[IX(i, j)] = 0;
        // Assume that densities below tolerance are too small to convey
        // any appreciable momentum to the dye
        m_vx[IX(i, j)] = 0;       
        m_vy[IX(i, j)] = 0; 
      }
      // else if the density and velocity of this cell is close enough to its neighbours',
      // don't simulate it
      else
      {
        abs_vx = std::abs(m_vx[IX(i, j)]);
        abs_vy = std::abs(m_vy[IX(i, j)]);
        if (
              std::abs(m_density[IX(i, j)] - m_density[IX(i + 1, j)]) < 1e-2 &&
              std::abs(m_density[IX(i, j)] - m_density[IX(i - 1, j)]) < 1e-2 &&
              std::abs(m_density[IX(i, j)] - m_density[IX(i, j + 1)]) < 1e-2 &&
              std::abs(m_density[IX(i, j)] - m_density[IX(i, j - 1)]) < 1e-2
              &&
              std::abs(abs_vx - std::abs(m_vx[IX(i + 1, j)]))         < 1e-3 &&
              std::abs(abs_vx - std::abs(m_vx[IX(i - 1, j)]))         < 1e-3 &&
              std::abs(abs_vx - std::abs(m_vx[IX(i, j + 1)]))         < 1e-3 &&
              std::abs(abs_vx - std::abs(m_vx[IX(i, j - 1)]))         < 1e-3 &&
              std::abs(abs_vy - std::abs(m_vy[IX(i + 1, j)]))         < 1e-3 &&
              std::abs(abs_vy - std::abs(m_vy[IX(i - 1, j)]))         < 1e-3 &&
              std::abs(abs_vy - std::abs(m_vy[IX(i, j + 1)]))         < 1e-3 &&
              std::abs(abs_vy - std::abs(m_vy[IX(i, j - 1)]))         < 1e-3
            )
        {
          m_surface[IX(i, j)] = 0;
        }
        else
        { 
          m_surface[IX(i, j    )] = 1;
          m_surface[IX(i, j + 1)] = 1;
          m_surface[IX(i, j - 1)] = 1;
          m_surface[IX(i + 1, j)] = 1;
          m_surface[IX(i - 1, j)] = 1; 
        }
      }
    }
  }
}

void FluidGrid::diffuse(uint8_t _b, float* _x, float* _x0, float _diff, float _dt)
{
  auto a = _dt * _diff * (Nx - 2) * (Ny - 2);
  linearSolve(_b, _x, _x0, a, (1.f + 6.f * a));
}

// prevents fluid "leaking" by treating all outer cells as walls
void FluidGrid::setBoundary(uint8_t _b, float* _x)
{
  int i;
  switch(_b)
  {
    case 0:
      for(i = 1; i < Ny - 1; i++) 
      {
        // Left
        _x[IX(0     , i)] = _x[IX(1     , i)];
        // Right
        _x[IX(Nx - 1, i)] = _x[IX(Nx - 2, i)];
      }
      for(i = 1; i < Nx - 1; i++) 
      {
        // Top
        _x[IX(i, 0    )] = _x[IX(i, 1    )];
        // Bottom
        _x[IX(i, Ny - 1)] = _x[IX(i, Ny - 2)];
      }
    case 1:
      for(i = 1; i < Ny - 1; i++) 
      {
        // Left
        _x[IX(0    , i)]  = -_x[IX(1     , i)];
        // Right
        _x[IX(Nx - 1, i)] = -_x[IX(Nx - 2, i)];
      }
      for(i = 1; i < Nx - 1; i++) 
      {
        // Top
        _x[IX(i, 0    )] = _x[IX(i, 1    )];
        // Bottom
        _x[IX(i, Ny - 1)] = _x[IX(i, Ny - 2)];
      }
    case 2:
      for(i = 1; i < Ny - 1; i++) 
      {
        // Left
        _x[IX(0     , i)] = _x[IX(1,      i)];
        // Right
        _x[IX(Nx - 1, i)] = _x[IX(Nx - 2, i)];
      }
      for(i = 1; i < Nx - 1; i++) 
      {
        // Top
        _x[IX(i, 0     )] = -_x[IX(i, 1     )];
        // Bottom
        _x[IX(i, Ny - 1)] = -_x[IX(i, Ny - 2)];
      }
  }
  
  _x[IX(0,      0     )] = .5f * (_x[IX(1,      0     )] + _x[IX(0,      1     )]);
  _x[IX(0,      Ny - 1)] = .5f * (_x[IX(1,      Ny - 1)] + _x[IX(0,      Ny - 2)]);
  _x[IX(Nx - 1, 0     )] = .5f * (_x[IX(Nx - 2, 0     )] + _x[IX(Nx - 1, 1     )]);
  _x[IX(Nx - 1, Ny - 1)] = .5f * (_x[IX(Nx - 2, Ny - 1)] + _x[IX(Nx - 1, Ny - 2)]);
}

void FluidGrid::linearSolve(uint8_t _b, float* _x, float* _x0, float _a, float _c)
{
  int i, j;
  auto cRecip = 1.f / _c;
  for (int k = 0; k < ITER; k++) 
  {
    FOR_EACH_CELL
      _x[IX(i, j)] = (_x0[IX(i, j)] + _a
                   * (_x[IX(i + 1, j    )]
                   +  _x[IX(i - 1, j    )]
                   +  _x[IX(i    , j + 1)]
                   +  _x[IX(i    , j - 1)])) * cRecip;
    END_FOR
    setBoundary(_b, _x);
  }
}

// ensures ingoing fluid equals outgoing fluid for every cell 
void FluidGrid::project(float* _velocX, float* _velocY, float* _p, float* _div)
{
  int i, j;
  auto NxRecip = 1.f / Nx;
  auto NyRecip = 1.f / Ny;
  FOR_EACH_CELL
    _div[IX(i, j)] = -0.5f*((_velocX[IX(i + 1, j    )] - _velocX[IX(i - 1, j    )]) * NxRecip
                          + (_velocY[IX(i    , j + 1)] - _velocY[IX(i    , j - 1)]) * NyRecip);
    _p[IX(i, j)] = 0;
  END_FOR
  setBoundary(0, _div); 
  setBoundary(0, _p);
  linearSolve(0, _p, _div, 1, 6);
  
  FOR_EACH_CELL
      _velocX[IX(i, j)] -= 0.5f * (_p[IX(i + 1, j)] - _p[IX(i - 1, j)]) * Nx;
      _velocY[IX(i, j)] -= 0.5f * (_p[IX(i, j + 1)] - _p[IX(i, j - 1)]) * Ny;
  END_FOR
  setBoundary(1, _velocX);
  setBoundary(2, _velocY);

}

void FluidGrid::advect(uint8_t _b, float* _d, float* _d0,  float* _velocX, float* _velocY, float _dt)
{  
  int i, j, i0, i1, j0, j1;
  float x, y, s0, s1, t0, t1;
  auto dtx = _dt * (Nx - 2);
  auto dty = _dt * (Ny - 2);

  FOR_EACH_CELL
    x = i - (dtx * _velocX[IX(i, j)]);
    y = j - (dty * _velocY[IX(i, j)]);
    x = std::max(.5f, std::min(x, Nx - 1.5f));
    y = std::max(.5f, std::min(y, Ny - 1.5f));
    i0 = int(x); 
    i1 = i0 + 1;
    j0 = int(y);
    j1 = j0 + 1;
    s1 = x - i0; 
    s0 = 1 - s1; 
    t1 = y - j0; 
    t0 = 1 - t1;
    _d[IX(i, j)] = s0 * (t0 * _d0[IX(i0, j0)] + t1 * _d0[IX(i0, j1)])
                 + s1 * (t0 * _d0[IX(i1, j0)] + t1 * _d0[IX(i1, j1)]);
  END_FOR
  setBoundary(_b, _d);
}

void FluidGrid::step(float _dt)
{
  // Optimisation stuff

  if (m_optimised)
    prune();

  // Velocity stuff

  // Diffuse the two velocity components
  // Pretty sure cell culling wouldn't effect this
  diffuse(1, m_vx0, m_vx, m_viscosity, _dt);
  diffuse(2, m_vy0, m_vy, m_viscosity, _dt);
  // Fix up velocities so they keep things incompressible
  project(m_vx0, m_vy0, m_vx, m_vy);
  // Move velocities around according to the velocities of the fluid
  advect(1, m_vx, m_vx0, m_vx0, m_vy0, _dt);
  advect(2, m_vy, m_vy0, m_vx0, m_vy0, _dt);
  // Fix up the velocities again
  project(m_vx, m_vy, m_vx0, m_vy0);

  // Density stuff

  // Diffuse the dye
  diffuse(0, m_s, m_density, m_diffusion, _dt);
  // Move the dye around according to the velocities
  advect(0, m_density, m_s, m_vx, m_vy, _dt);
}


int FluidGrid::getDensity(int _i)
{
  return std::min(int(m_density[_i]), 255);
}

int FluidGrid::isOptimised()
{
  return m_optimised;
}

int FluidGrid::getSurfaceFlag(int _x, int _y)
{
  return m_surface[IX(_x, _y)];
}

int FluidGrid::getSurfaceFlag(int _i)
{
  return m_surface[_i];
}

int FluidGrid::getLeftEdge()
{
  return m_leftEdge;
}

int FluidGrid::getRightEdge()
{
  return m_rightEdge;
}

int FluidGrid::getTopEdge()
{
  return m_topEdge;
}

int FluidGrid::getBottomEdge()
{
  return m_bottomEdge;
}






#pragma once

#include <string>

#include "mesh_structured.h"
#include "physics_flowparameters.h"

// ---------------------------------------------------------------------------
// Test cases
//
//  applyTestCase() tags the boundary faces, writes the initial condition
//  into mesh.cells[].U, finalises fp (buildDerived) and returns the end
//  time of the case.  Every case except Channel overwrites fp with its own
//  inflow / reference state.
//
//    Channel    Inlet left, outlet right, slip walls.  Starts at rest with
//               fp's density and pressure.  Runs until stopped.
//    Sod        Sod (1978) shock tube, membrane at Lx/2.  t_end = 0.2
//    ShuOsher   Shu & Osher (1989): Mach-3 shock running into a sinusoidal
//               density field.  Canonical domain x in [-5, 5], mapped to
//               [0, 10] (use Lx = 10).  t_end = 1.8
//    RiemannN   Lax & Liu (1998) 2-D Riemann problem, configuration N.
//               Four constant quadrants meeting at the domain centre,
//               transmissive (zero-gradient) boundaries on all sides.
//               Canonical domain: unit square.
//    SmoothWave Gaussian density bump advected at u = 2 with p = 1
//               (supersonic, so the inlet and outlet are exact).  Smooth
//               exact solution smoothWaveDensity(x, t, Lx): the order-of-
//               accuracy test.  t_end = 0.2 Lx
// ---------------------------------------------------------------------------

enum class TestCase
{
    Channel,
    Sod,
    ShuOsher,
    Riemann3,    // four shocks                       t_end = 0.3
    Riemann4,    // four shocks                       t_end = 0.25
    Riemann6,    // four contacts (vortex sheets)     t_end = 0.3
    Riemann12,   // two shocks + two contacts         t_end = 0.25
    SmoothWave
};

// Exact density of SmoothWave: 1 + 0.2 exp(-((x - 2t - 0.3 Lx) / (0.08 Lx))^2)
double smoothWaveDensity(double x, double t, double Lx);

// Domain length Shu-Osher needs: its shock and density wave sit at absolute
// positions of the canonical domain x in [-5, 5]
constexpr double SHU_OSHER_LX = 10.0;

// Returns t_end (1e30 for Channel, which has no natural end time)
double applyTestCase(TestCase tc, StructuredMesh& mesh, FlowParameters& fp,
                     bool subsonic);

// Short identifier used in file metadata and on the command line:
// channel, sod, shu-osher, riemann3, riemann4, riemann6, riemann12
const char* testCaseId(TestCase tc);

// Inverse of testCaseId(); returns false for an unknown identifier
bool testCaseFromId(const std::string& id, TestCase& tc);

// Default mesh of each case (the canonical domain where there is one)
struct MeshSize { int Nx; int Ny; double Lx; double Ly; };
MeshSize canonicalMesh(TestCase tc);

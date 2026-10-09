#include "physics_testcases.h"
#include "physics_eulerphysics.h"

#include <cmath>

namespace {

void setCell(Cell& c, const PrimitiveState& w)
{
    c.U = EulerPhysics::primitiveToConserved(w.rho, w.u, w.v, w.p);
}

void setFlow(FlowParameters& fp, const PrimitiveState& w)
{
    fp.rho_inf = w.rho;
    fp.u_inf   = w.u;
    fp.v_inf   = w.v;
    fp.p_inf   = w.p;
}

// Lax & Liu (1998) initial data {rho, u, v, p} per quadrant:
// q[0] upper right, q[1] upper left, q[2] lower left, q[3] lower right.
struct RiemannConfig { PrimitiveState q[4]; double tEnd; };

RiemannConfig riemannConfig(TestCase tc)
{
    switch(tc)
    {
    case TestCase::Riemann3:
        return {{{1.5,    0.0,    0.0,    1.5  },
                 {0.5323, 1.206,  0.0,    0.3  },
                 {0.138,  1.206,  1.206,  0.029},
                 {0.5323, 0.0,    1.206,  0.3  }}, 0.3};
    case TestCase::Riemann4:
        return {{{1.1,    0.0,    0.0,    1.1  },
                 {0.5065, 0.8939, 0.0,    0.35 },
                 {1.1,    0.8939, 0.8939, 1.1  },
                 {0.5065, 0.0,    0.8939, 0.35 }}, 0.25};
    case TestCase::Riemann6:
        return {{{1.0,    0.75,  -0.5,    1.0  },
                 {2.0,    0.75,   0.5,    1.0  },
                 {1.0,   -0.75,   0.5,    1.0  },
                 {3.0,   -0.75,  -0.5,    1.0  }}, 0.3};
    default: // Riemann12
        return {{{0.5313, 0.0,    0.0,    0.4  },
                 {1.0,    0.7276, 0.0,    1.0  },
                 {0.8,    0.0,    0.0,    1.0  },
                 {1.0,    0.0,    0.7276, 1.0  }}, 0.25};
    }
}

} // anonymous namespace

double applyTestCase(TestCase tc, StructuredMesh& mesh, FlowParameters& fp,
                     bool subsonic)
{
    double tEnd = 1e30;

    switch(tc)
    {
    case TestCase::Channel:
        mesh.setBoundaries(subsonic ? BoundaryType::InletSubsonic  : BoundaryType::Inlet,
                           subsonic ? BoundaryType::OutletSubsonic : BoundaryType::Outlet,
                           BoundaryType::Wall);
        // Uniform state at rest; the inlet drives the flow
        for(Cell& c : mesh.cells)
            setCell(c, {fp.rho_inf, 0.0, 0.0, fp.p_inf});
        break;

    case TestCase::Sod:
    {
        const PrimitiveState L{1.0, 0.0, 0.0, 1.0}, R{0.125, 0.0, 0.0, 0.1};
        mesh.setBoundaries(BoundaryType::ZeroGradient, BoundaryType::ZeroGradient,
                           BoundaryType::Wall);
        for(Cell& c : mesh.cells)
            setCell(c, c.cx < 0.5 * mesh.Lx ? L : R);
        setFlow(fp, L);
        fp.p_back = R.p;
        tEnd = 0.2;
        break;
    }

    case TestCase::ShuOsher:
    {
        // Post-shock state of a Mach-3 shock into (rho, u, p) = (1, 0, 1).
        // It enters through the supersonic inlet; the outlet is never reached.
        const PrimitiveState post{3.857143, 2.629369, 0.0, 10.33333};
        mesh.setBoundaries(BoundaryType::Inlet, BoundaryType::Outlet,
                           BoundaryType::Wall);
        for(Cell& c : mesh.cells)
        {
            const double x = c.cx - 0.5 * SHU_OSHER_LX;
            setCell(c, x < -4.0 ? post
                                : PrimitiveState{1.0 + 0.2 * std::sin(5.0 * x), 0.0, 0.0, 1.0});
        }
        setFlow(fp, post);
        tEnd = 1.8;
        break;
    }

    case TestCase::SmoothWave:
    {
        // The bump is ~1e-16 at both ends, so the constant inlet state and
        // the zero-gradient outlet are exact
        mesh.setBoundaries(BoundaryType::Inlet, BoundaryType::Outlet,
                           BoundaryType::Wall);
        for(Cell& c : mesh.cells)
            setCell(c, {smoothWaveDensity(c.cx, 0.0, mesh.Lx), 2.0, 0.0, 1.0});
        setFlow(fp, {1.0, 2.0, 0.0, 1.0});
        tEnd = 0.2 * mesh.Lx;
        break;
    }

    default: // 2-D Riemann problems
    {
        const RiemannConfig cfg = riemannConfig(tc);
        mesh.setBoundaries(BoundaryType::ZeroGradient, BoundaryType::ZeroGradient,
                           BoundaryType::ZeroGradient);
        for(Cell& c : mesh.cells)
        {
            const bool right = c.cx >= 0.5 * mesh.Lx;
            const bool top   = c.cy >= 0.5 * mesh.Ly;
            setCell(c, cfg.q[top ? (right ? 0 : 1) : (right ? 3 : 2)]);
        }
        setFlow(fp, cfg.q[0]);
        tEnd = cfg.tEnd;
        break;
    }
    }

    fp.buildDerived();
    return tEnd;
}

double smoothWaveDensity(double x, double t, double Lx)
{
    const double s = (x - 2.0 * t - 0.3 * Lx) / (0.08 * Lx);
    return 1.0 + 0.2 * std::exp(-s * s);
}

namespace {

struct CaseEntry { TestCase tc; const char* id; MeshSize mesh; };

const CaseEntry CASES[] = {
    {TestCase::Channel,   "channel",   {100,  20, 1.0,          0.2}},
    {TestCase::Sod,       "sod",       {200,   1, 1.0,          1.0}},
    {TestCase::ShuOsher,  "shu-osher", {400,   1, SHU_OSHER_LX, 1.0}},
    {TestCase::Riemann3,  "riemann3",  {200, 200, 1.0,          1.0}},
    {TestCase::Riemann4,  "riemann4",  {200, 200, 1.0,          1.0}},
    {TestCase::Riemann6,  "riemann6",  {200, 200, 1.0,          1.0}},
    {TestCase::Riemann12, "riemann12", {200, 200, 1.0,          1.0}},
    {TestCase::SmoothWave, "smooth-wave", {200, 1, 1.0,         1.0}},
};

const CaseEntry& entry(TestCase tc)
{
    for(const CaseEntry& e : CASES)
        if(e.tc == tc) return e;
    return CASES[0];   // unreachable: every enum value is listed
}

} // anonymous namespace

const char* testCaseId(TestCase tc)
{
    return entry(tc).id;
}

bool testCaseFromId(const std::string& id, TestCase& tc)
{
    for(const CaseEntry& e : CASES)
        if(id == e.id) { tc = e.tc; return true; }
    return false;
}

MeshSize canonicalMesh(TestCase tc)
{
    return entry(tc).mesh;
}

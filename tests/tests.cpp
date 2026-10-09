// Headless regression checks (no Qt). Build the `solver_tests` target and
// run it directly or through `ctest`. Exit code 0 = all checks passed.

#include "physics_solver.h"
#include "physics_testcases.h"
#include "sod_exact.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>

static int failures = 0;

static void check(bool ok, const char* what)
{
    std::printf("%s  %s\n", ok ? "PASS" : "FAIL", what);
    if(!ok) ++failures;
}

// The 8 solver configurations: bit 0 = HLLC, bit 1 = MUSCL, bit 2 = SSP-RK2
static void configure(EulerSolver& solver, int cfg)
{
    solver.setFluxScheme(cfg & 1 ? EulerSolver::FluxScheme::HLLC
                                 : EulerSolver::FluxScheme::Rusanov);
    solver.setReconstructionScheme(cfg & 2 ? EulerSolver::ReconstructionScheme::MUSCL
                                           : EulerSolver::ReconstructionScheme::PiecewiseConstant);
    solver.setTimeScheme(cfg & 4 ? EulerSolver::TimeScheme::RK2
                                 : EulerSolver::TimeScheme::ForwardEuler);
}

static double maxAbs(const EulerState& U)
{
    return std::max({std::abs(U.rho), std::abs(U.rho_u), std::abs(U.rho_v), std::abs(U.rho_E)});
}

// Domain totals of the conserved variables, sum of U * area
static EulerState totals(const StructuredMesh& mesh)
{
    EulerState sum;
    for(const Cell& c : mesh.cells)
        sum += c.U * c.area;
    return sum;
}

// Runs a test case to its end time with HLLC + MUSCL + SSP-RK2
static void runCase(TestCase tc, StructuredMesh& mesh)
{
    FlowParameters fp;
    const double tEnd = applyTestCase(tc, mesh, fp, false);
    EulerSolver solver(&mesh, fp);
    configure(solver, 7);
    solver.setTEnd(tEnd);
    while(solver.simTime() < solver.tEnd())
        solver.step();
}

static void runTests()
{
    // 1. Exact Riemann solver reproduces the textbook Sod star state (Toro)
    SodExact sod({1.0, 0.0, 1.0}, {0.125, 0.0, 0.1});
    check(std::abs(sod.pStar() - 0.30313) < 1e-4 &&
          std::abs(sod.uStar() - 0.92745) < 1e-4, "Sod p* and u*");

    // 2. Mirrored Sod (rarefaction on the RIGHT) is the mirror image of Sod
    SodExact mirror({0.125, 0.0, 0.1}, {1.0, 0.0, 1.0});
    bool symmetric = true;
    for(double x : {0.1, 0.3, 0.45, 0.55, 0.7, 0.9})
    {
        SodExact::Point a = sod.sample(x, 0.2, 0.5);
        SodExact::Point b = mirror.sample(1.0 - x, 0.2, 0.5);
        symmetric = symmetric && std::abs(a.rho - b.rho) < 1e-9
                              && std::abs(a.u   + b.u)   < 1e-9
                              && std::abs(a.p   - b.p)   < 1e-9;
    }
    check(symmetric, "right rarefaction fan mirrors the left one");

    // 3. Switching cases on one mesh retags every side (walls included)
    {
        StructuredMesh mesh(8, 8, 1.0, 1.0);
        FlowParameters fp;
        applyTestCase(TestCase::Riemann3, mesh, fp, false);   // all zero-gradient
        applyTestCase(TestCase::Channel,  mesh, fp, false);
        bool ok = true;
        for(const Face& f : mesh.faces)
        {
            if(f.boundaryType == BoundaryType::Internal) continue;
            const BoundaryType expected = f.nx < 0.0 ? BoundaryType::Inlet
                                        : f.nx > 0.0 ? BoundaryType::Outlet
                                                     : BoundaryType::Wall;
            ok = ok && f.boundaryType == expected;
        }
        check(ok, "boundary re-tagging between test cases");
    }

    // 4. Sod: mass conservation and L1 density error vs the exact solution
    {
        StructuredMesh mesh(100, 1, 1.0, 0.1);
        runCase(TestCase::Sod, mesh);
        double mass = 0.0, l1 = 0.0;
        for(const Cell& c : mesh.cells)
        {
            mass += c.U.rho * c.area;
            l1   += std::abs(c.U.rho - sod.sample(c.cx, 0.2, 0.5).rho);
        }
        l1 /= mesh.cells.size();
        std::printf("      Sod N=100: L1(rho) = %.6f\n", l1);
        check(std::abs(mass - 0.5625 * 0.1) < 1e-12, "Sod mass conservation");
        check(l1 < 0.01, "Sod L1 density error");
    }

    // 5. Shu-Osher: shock position (x ~ 2.40 on converged grids) and the
    //    density wave ahead of it left untouched (HLLC keeps stationary
    //    contacts exact; Rusanov would smear it)
    {
        StructuredMesh mesh(400, 1, 10.0, 1.0);
        runCase(TestCase::ShuOsher, mesh);
        double xShock = -5.0;
        for(const Cell& c : mesh.cells)
            if(EulerPhysics::pressure(c.U) > 1.5)
                xShock = std::max(xShock, c.cx - 5.0);
        double drift = 0.0;
        for(const Cell& c : mesh.cells)
        {
            const double x = c.cx - 5.0;
            if(x > xShock + 0.5)
                drift = std::max(drift, std::abs(c.U.rho - (1.0 + 0.2 * std::sin(5.0 * x))));
        }
        std::printf("      Shu-Osher N=400: shock at x = %.4f\n", xShock);
        check(std::abs(xShock - 2.40) < 0.05, "Shu-Osher shock position");
        check(drift < 1e-12, "Shu-Osher undisturbed region ahead of the shock");
    }

    // 6. 2-D Riemann config 3 is symmetric about the diagonal y = x:
    //    rho(i,j) = rho(j,i) and rho*u(i,j) = rho*v(j,i)
    {
        const int N = 64;
        StructuredMesh mesh(N, N, 1.0, 1.0);
        runCase(TestCase::Riemann3, mesh);
        double asym = 0.0;
        for(int j = 0; j < N; ++j)
            for(int i = 0; i < N; ++i)
            {
                const EulerState& a = mesh.cells[mesh.cellIndex(i, j)].U;
                const EulerState& b = mesh.cells[mesh.cellIndex(j, i)].U;
                asym = std::max({asym, std::abs(a.rho   - b.rho),
                                       std::abs(a.rho_u - b.rho_v),
                                       std::abs(a.rho_E - b.rho_E)});
            }
        check(asym < 1e-10, "2-D Riemann config 3 diagonal symmetry");
    }

    // 7. Every 2-D Riemann configuration reaches t_end with a physical state
    {
        bool ok = true;
        for(TestCase tc : {TestCase::Riemann3, TestCase::Riemann4,
                           TestCase::Riemann6, TestCase::Riemann12})
        {
            StructuredMesh mesh(32, 32, 1.0, 1.0);
            runCase(tc, mesh);
            for(const Cell& c : mesh.cells)
                ok = ok && c.U.rho > 0.0 && EulerPhysics::pressure(c.U) > 0.0;
        }
        check(ok, "2-D Riemann configs 3/4/6/12 stay positive");
    }

    // 8. HLLC stays finite when a state has p <= 0 (possible in a MUSCL
    //    face state); without the pressure floor p = 0 gave NaN
    {
        const EulerState UL = EulerPhysics::primitiveToConserved(1.0, 0.0, 0.0, 1.0);
        bool finite = true;
        for(double pR : {0.0, -0.5})
        {
            const EulerState UR = EulerPhysics::primitiveToConserved(1.0, 0.0, 0.0, pR);
            const EulerState F  = EulerPhysics::hllcFlux(UL, UR, 1.0, 0.0);
            finite = finite && std::isfinite(F.rho)   && std::isfinite(F.rho_u)
                            && std::isfinite(F.rho_v) && std::isfinite(F.rho_E);
        }
        check(finite, "HLLC flux finite for p <= 0");
    }

    // 9. An unstable run (CFL twice the limit) stops with an error instead
    //    of carrying on with a non-physical field
    {
        StructuredMesh mesh(100, 1, 1.0, 0.1);
        FlowParameters fp;
        const double tEnd = applyTestCase(TestCase::Sod, mesh, fp, false);
        EulerSolver solver(&mesh, fp);
        solver.setCFL(4.0);
        solver.setTEnd(tEnd);
        bool stopped = false;
        try
        {
            while(solver.simTime() < solver.tEnd())
                solver.step();
        }
        catch(const std::runtime_error& e)
        {
            std::printf("      %s\n", e.what());
            stopped = true;
        }
        check(stopped, "blow-up stops the run");
    }

    // 10. Both fluxes are consistent, F(U,U,n) = F(U)·n, and antisymmetric,
    //     F(UL,UR,n) = -F(UR,UL,-n), for any face orientation
    {
        const EulerState UL = EulerPhysics::primitiveToConserved(1.0,    0.3, -0.2, 1.0);
        const EulerState UR = EulerPhysics::primitiveToConserved(0.125, -0.4,  0.5, 0.1);
        const double s = std::sqrt(0.5);
        const double normals[][2] = {{1, 0}, {0, 1}, {-1, 0}, {s, s}, {0.6, -0.8}};
        double consistency = 0.0, antisymmetry = 0.0;
        for(bool hllc : {false, true})
        {
            auto flux = [hllc](const EulerState& a, const EulerState& b, double nx, double ny)
            {
                return hllc ? EulerPhysics::hllcFlux(a, b, nx, ny)
                            : EulerPhysics::rusanovFlux(a, b, nx, ny);
            };
            for(const auto& n : normals)
            {
                consistency  = std::max(consistency,
                    maxAbs(flux(UL, UL, n[0], n[1]) - EulerPhysics::physicalFluxNormal(UL, n[0], n[1])));
                antisymmetry = std::max(antisymmetry,
                    maxAbs(flux(UL, UR, n[0], n[1]) + flux(UR, UL, -n[0], -n[1])));
            }
        }
        check(consistency  < 1e-13, "fluxes consistent: F(U,U) = F(U).n");
        check(antisymmetry < 1e-13, "fluxes antisymmetric under n -> -n");
    }

    // 11. A uniform oblique stream with far-field boundaries stays uniform in
    //     all 8 configurations (every cell is closed, boundary normals point
    //     outward, x- and y-faces are assembled consistently)
    {
        double drift = 0.0;
        for(int cfg = 0; cfg < 8; ++cfg)
        {
            StructuredMesh mesh(12, 8, 1.5, 1.0);
            mesh.setBoundaries(BoundaryType::Farfield, BoundaryType::Farfield,
                               BoundaryType::Farfield);
            FlowParameters fp;
            fp.rho_inf = 1.0;  fp.u_inf = 0.5;  fp.v_inf = -0.3;  fp.p_inf = 1.0;
            fp.buildDerived();
            for(Cell& c : mesh.cells)
                c.U = fp.U_inlet;
            EulerSolver solver(&mesh, fp);
            configure(solver, cfg);
            for(int n = 0; n < 20; ++n)
                solver.step();
            for(const Cell& c : mesh.cells)
                drift = std::max(drift, maxAbs(c.U - fp.U_inlet));
        }
        check(drift < 1e-13, "uniform stream preserved by all 8 configurations");
    }

    // 12. The Channel boundaries (supersonic or characteristic subsonic inlet
    //     and outlet, slip walls) hold a uniform stream when p_back = p
    {
        double drift = 0.0;
        for(double u : {0.5, 2.0})   // subsonic and supersonic (c = 1.18)
            for(int cfg = 0; cfg < 8; ++cfg)
            {
                StructuredMesh mesh(12, 4, 1.0, 0.5);
                FlowParameters fp;
                fp.rho_inf = 1.0;  fp.u_inf = u;  fp.v_inf = 0.0;
                fp.p_inf   = 1.0;  fp.p_back = 1.0;
                applyTestCase(TestCase::Channel, mesh, fp, u < 1.0);
                for(Cell& c : mesh.cells)   // start from the stream, not at rest
                    c.U = fp.U_inlet;
                EulerSolver solver(&mesh, fp);
                configure(solver, cfg);
                for(int n = 0; n < 20; ++n)
                    solver.step();
                for(const Cell& c : mesh.cells)
                    drift = std::max(drift, maxAbs(c.U - fp.U_inlet));
            }
        check(drift < 1e-12, "channel boundaries hold a uniform stream");
    }

    // 13. In a closed box (slip walls on every side) mass and energy are
    //     conserved to round-off while the 2-D Riemann waves reflect
    {
        double err = 0.0;
        for(int cfg : {0, 7})   // Rusanov PC FE and HLLC MUSCL SSP-RK2
        {
            StructuredMesh mesh(32, 32, 1.0, 1.0);
            FlowParameters fp;
            const double tEnd = applyTestCase(TestCase::Riemann3, mesh, fp, false);
            mesh.setBoundaries(BoundaryType::Wall, BoundaryType::Wall, BoundaryType::Wall);
            const EulerState before = totals(mesh);
            EulerSolver solver(&mesh, fp);
            configure(solver, cfg);
            solver.setTEnd(2.0 * tEnd);   // long enough to hit the walls
            while(solver.simTime() < solver.tEnd())
                solver.step();
            const EulerState after = totals(mesh);
            err = std::max({err, std::abs(after.rho   - before.rho)   / before.rho,
                                 std::abs(after.rho_E - before.rho_E) / before.rho_E});
        }
        check(err < 1e-12, "closed box conserves mass and energy");
    }

    // 14. Order of accuracy on a smooth density wave advected at u = 2.
    //     The flow is supersonic, so nothing travels upstream from the
    //     outlet; the error is measured away from the inlet, whose constant
    //     state does not follow the wave.
    {
        auto l1Error = [](int N, int cfg)
        {
            const double pi = std::acos(-1.0);
            StructuredMesh mesh(N, 1, 1.0, 1.0);
            mesh.setBoundaries(BoundaryType::Inlet, BoundaryType::Outlet, BoundaryType::Wall);
            FlowParameters fp;
            fp.rho_inf = 1.0;  fp.u_inf = 2.0;  fp.v_inf = 0.0;  fp.p_inf = 1.0;
            fp.buildDerived();
            for(Cell& c : mesh.cells)
                c.U = EulerPhysics::primitiveToConserved(
                    1.0 + 0.2 * std::sin(2.0 * pi * c.cx), 2.0, 0.0, 1.0);
            EulerSolver solver(&mesh, fp);
            configure(solver, cfg);
            solver.setTEnd(0.1);
            while(solver.simTime() < solver.tEnd())
                solver.step();
            double l1 = 0.0;
            int    n  = 0;
            for(const Cell& c : mesh.cells)
                if(c.cx > 0.4 && c.cx < 0.9)
                {
                    l1 += std::abs(c.U.rho - (1.0 + 0.2 * std::sin(2.0 * pi * (c.cx - 0.2))));
                    ++n;
                }
            return l1 / n;
        };
        const double orderPC    = std::log2(l1Error(200, 5) / l1Error(400, 5));   // HLLC PC RK2
        const double orderMUSCL = std::log2(l1Error(200, 7) / l1Error(400, 7));   // HLLC MUSCL RK2
        std::printf("      observed order: PC %.2f, MUSCL %.2f\n", orderPC, orderMUSCL);
        check(orderPC > 0.9 && orderPC < 1.1, "piecewise constant is first order");
        check(orderMUSCL > 1.8, "MUSCL + SSP-RK2 is second order on smooth flow");
    }
}

int main()
{
    // A run that blows up throws; report it instead of aborting
    try
    {
        runTests();
    }
    catch(const std::exception& e)
    {
        std::printf("FAIL  unexpected exception: %s\n", e.what());
        ++failures;
    }
    return failures ? 1 : 0;
}

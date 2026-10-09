// ---------------------------------------------------------------------------
// euler_cli — headless runner
//
//  Runs one test case with the chosen schemes and writes CSV files in the
//  same format as the GUI export, plus run.json.  Meant for scripted studies
//  (see analysis/).  Run with --help for the options.
// ---------------------------------------------------------------------------

#include "io_csv.h"
#include "physics_solver.h"
#include "physics_testcases.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

const char* USAGE = R"(Usage: euler_cli --case NAME --out DIR [options]

Runs one test case and writes CSV files (same format as the GUI export)
and run.json into DIR.

  --case NAME          channel | sod | shu-osher | riemann3 | riemann4 |
                       riemann6 | riemann12 | smooth-wave
  --out DIR            output directory (created if needed)

  --nx N  --ny N       cells                     [canonical mesh of the case]
  --lx L  --ly L       domain size               [canonical mesh of the case]
  --flux F             rusanov | hllc            [hllc]
  --recon R            pc | muscl                [muscl]
  --time T             fe | rk2                  [rk2]
  --cfl C              CFL number, about twice the classical Courant number [0.5]
  --t-end T            end time                  [the case's own]
  --max-iter N         step limit                [100000]
  --snapshot-every K   also write snapshot_iter_*.csv every K steps [0 = off]
  --probe-j J          row of the line probe     [Ny / 2]

  Channel only (the other cases set their own flow state):
  --rho R  --u U  --v V  --p P  --p-back PB
                       free stream and back pressure [1.225 0 0 101325 101325]
  --subsonic           characteristic subsonic inlet and outlet

Output: final.csv (all cells at the end), probe_jNN.csv (row J at the end),
residuals.csv, run.json and the periodic snapshots if requested.
Exit code: 0 = finished, 1 = invalid arguments or I/O error,
2 = the solution became non-physical (blow-up).
)";

const std::set<std::string> VALUE_OPTIONS = {
    "--case", "--out", "--nx", "--ny", "--lx", "--ly", "--flux", "--recon", "--time",
    "--cfl", "--t-end", "--max-iter", "--snapshot-every", "--probe-j",
    "--rho", "--u", "--v", "--p", "--p-back"};

const std::set<std::string> FLAG_OPTIONS = {"--subsonic", "--help"};

// ---------------------------------------------------------------------------
//  Arguments: "--key value" pairs and value-less flags
// ---------------------------------------------------------------------------
using Args = std::map<std::string, std::string>;

Args parseArgs(int argc, char** argv)
{
    Args args;
    for(int i = 1; i < argc; ++i)
    {
        const std::string key = argv[i];
        if(FLAG_OPTIONS.count(key))
            args[key] = "1";
        else if(!VALUE_OPTIONS.count(key))
            throw std::invalid_argument("unknown option '" + key + "'");
        else if(i + 1 >= argc)
            throw std::invalid_argument("missing value for " + key);
        else
            args[key] = argv[++i];
    }
    return args;
}

double number(const Args& args, const std::string& key, double fallback)
{
    const auto it = args.find(key);
    if(it == args.end()) return fallback;

    std::size_t used = 0;
    double value = 0.0;
    try { value = std::stod(it->second, &used); }
    catch(const std::exception&) { used = 0; }
    if(used == 0 || used != it->second.size() || !std::isfinite(value))
        throw std::invalid_argument("invalid number for " + key + ": '" + it->second + "'");
    return value;
}

int integer(const Args& args, const std::string& key, int fallback)
{
    const double value = number(args, key, fallback);
    if(value != static_cast<int>(value))
        throw std::invalid_argument(key + " must be an integer");
    return static_cast<int>(value);
}

std::string choice(const Args& args, const std::string& key, const std::string& fallback,
                   const std::set<std::string>& allowed)
{
    const auto it = args.find(key);
    const std::string value = (it == args.end()) ? fallback : it->second;
    if(!allowed.count(value))
        throw std::invalid_argument("invalid value for " + key + ": '" + value + "'");
    return value;
}

void require(bool ok, const std::string& what)
{
    if(!ok) throw std::invalid_argument(what);
}

// ---------------------------------------------------------------------------
//  Output helpers
// ---------------------------------------------------------------------------

// Writes one file; throws std::runtime_error if it cannot be written
template<class Write>
void writeFile(const fs::path& path, Write write)
{
    std::ofstream out(path);
    if(out)
        write(out);
    out.close();
    if(!out)
        throw std::runtime_error("cannot write " + path.string());
}

std::vector<EulerState> states(const StructuredMesh& mesh)
{
    std::vector<EulerState> U;
    U.reserve(mesh.cells.size());
    for(const Cell& c : mesh.cells)
        U.push_back(c.U);
    return U;
}

std::string jsonEscape(const std::string& s)
{
    std::string out;
    for(char c : s)
    {
        if(c == '"' || c == '\\') out += '\\';
        if(static_cast<unsigned char>(c) >= 0x20) out += c;
    }
    return out;
}

void writeRunJson(std::ostream& out, const RunInfo& r, int steps, double simTime,
                  double wallSeconds, const std::string& status, const std::string& message)
{
    out << std::setprecision(10)
        << "{\n"
        << "  \"case\": \""           << r.caseName       << "\",\n"
        << "  \"nx\": "               << r.Nx             << ",\n"
        << "  \"ny\": "               << r.Ny             << ",\n"
        << "  \"lx\": "               << r.Lx             << ",\n"
        << "  \"ly\": "               << r.Ly             << ",\n"
        << "  \"flux\": \""           << r.flux           << "\",\n"
        << "  \"reconstruction\": \"" << r.reconstruction << "\",\n"
        << "  \"time_scheme\": \""    << r.timeScheme     << "\",\n"
        << "  \"cfl\": "              << r.CFL            << ",\n"
        << "  \"t_end\": "            << r.tEnd           << ",\n"
        << "  \"max_iter\": "         << r.maxIter        << ",\n"
        << "  \"subsonic\": "         << (r.subsonic ? "true" : "false") << ",\n"
        << "  \"rho_inf\": "          << r.rho_inf        << ",\n"
        << "  \"u_inf\": "            << r.u_inf          << ",\n"
        << "  \"v_inf\": "            << r.v_inf          << ",\n"
        << "  \"p_inf\": "            << r.p_inf          << ",\n"
        << "  \"p_back\": "           << r.p_back         << ",\n"
        << "  \"steps\": "            << steps            << ",\n"
        << "  \"sim_time\": "         << simTime          << ",\n"
        << "  \"wall_time_s\": "      << wallSeconds      << ",\n"
        << "  \"status\": \""         << status           << "\",\n"
        << "  \"message\": \""        << jsonEscape(message) << "\"\n"
        << "}\n";
}

std::string snapshotName(int iter)
{
    char name[64];
    std::snprintf(name, sizeof name, "snapshot_iter_%07d.csv", iter);
    return name;
}

// ---------------------------------------------------------------------------
//  One run: returns the exit code
// ---------------------------------------------------------------------------
int run(const Args& args)
{
    // ---- Case and mesh ----------------------------------------------------
    require(args.count("--case"), "--case is required");
    require(args.count("--out"),  "--out is required");

    TestCase tc;
    require(testCaseFromId(args.at("--case"), tc),
            "unknown case '" + args.at("--case") + "'");

    const MeshSize m  = canonicalMesh(tc);
    const int    nx   = integer(args, "--nx", m.Nx);
    const int    ny   = integer(args, "--ny", m.Ny);
    const double lx   = number (args, "--lx", m.Lx);
    const double ly   = number (args, "--ly", m.Ly);
    require(nx >= 1 && ny >= 1, "--nx and --ny must be at least 1");
    require(lx > 0.0 && ly > 0.0, "--lx and --ly must be positive");

    // ---- Schemes and run control -----------------------------------------
    const std::string flux  = choice(args, "--flux",  "hllc",  {"rusanov", "hllc"});
    const std::string recon = choice(args, "--recon", "muscl", {"pc", "muscl"});
    const std::string time  = choice(args, "--time",  "rk2",   {"fe", "rk2"});
    const double cfl           = number (args, "--cfl", 0.5);
    const int    maxIter       = integer(args, "--max-iter", 100000);
    const int    snapshotEvery = integer(args, "--snapshot-every", 0);
    const int    probeJ        = integer(args, "--probe-j", ny / 2);
    const bool   subsonic      = args.count("--subsonic") > 0;
    require(cfl > 0.0, "--cfl must be positive");
    require(maxIter >= 1, "--max-iter must be at least 1");
    require(snapshotEvery >= 0, "--snapshot-every must not be negative");
    require(probeJ >= 0 && probeJ < ny, "--probe-j must be in [0, ny)");

    // ---- Flow state (Channel; the other cases overwrite it) ---------------
    FlowParameters fp;
    fp.rho_inf = number(args, "--rho",    fp.rho_inf);
    fp.u_inf   = number(args, "--u",      fp.u_inf);
    fp.v_inf   = number(args, "--v",      fp.v_inf);
    fp.p_inf   = number(args, "--p",      fp.p_inf);
    fp.p_back  = number(args, "--p-back", fp.p_back);
    require(fp.rho_inf > 0.0 && fp.p_inf > 0.0 && fp.p_back > 0.0,
            "--rho, --p and --p-back must be positive");

    // ---- Set up -------------------------------------------------------------
    StructuredMesh mesh(nx, ny, lx, ly);
    double tEnd = applyTestCase(tc, mesh, fp, subsonic);
    tEnd = number(args, "--t-end", tEnd);
    require(tEnd > 0.0, "--t-end must be positive");

    EulerSolver solver(&mesh, fp);
    solver.setFluxScheme(flux == "hllc" ? EulerSolver::FluxScheme::HLLC
                                        : EulerSolver::FluxScheme::Rusanov);
    solver.setReconstructionScheme(recon == "muscl"
                                       ? EulerSolver::ReconstructionScheme::MUSCL
                                       : EulerSolver::ReconstructionScheme::PiecewiseConstant);
    solver.setTimeScheme(time == "rk2" ? EulerSolver::TimeScheme::RK2
                                       : EulerSolver::TimeScheme::ForwardEuler);
    solver.setCFL(cfl);
    solver.setTEnd(tEnd);

    RunInfo info      = makeRunInfo(tc, subsonic, mesh, fp, solver);
    info.maxIter       = maxIter;
    info.snapshotEvery = snapshotEvery;

    const fs::path out = args.at("--out");
    fs::create_directories(out);

    auto writeSnapshot = [&](int iter) {
        writeFile(out / snapshotName(iter), [&](std::ostream& f) {
            writeSnapshotCsv(f, info, mesh, states(mesh), iter, solver.simTime());
        });
    };

    // ---- Time loop ----------------------------------------------------------
    std::string status = "max_iter", message;
    int steps = 0;
    if(snapshotEvery > 0)
        writeSnapshot(0);

    const auto start = std::chrono::steady_clock::now();
    while(steps < maxIter)
    {
        if(solver.simTime() >= solver.tEnd()) { status = "t_end"; break; }

        try
        {
            solver.step();
        }
        catch(const std::runtime_error& e)
        {
            status  = "blowup";
            message = e.what();
            break;
        }
        ++steps;

        if(snapshotEvery > 0 && steps % snapshotEvery == 0)
            writeSnapshot(steps);
    }
    if(status == "max_iter" && solver.simTime() >= solver.tEnd())
        status = "t_end";
    const double wall = std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - start).count();

    // ---- Output (also after a blow-up, to show where it started) ----------
    const std::vector<EulerState> U = states(mesh);
    writeFile(out / "final.csv", [&](std::ostream& f) {
        writeSnapshotCsv(f, info, mesh, U, steps, solver.simTime());
    });
    char probeName[32];
    std::snprintf(probeName, sizeof probeName, "probe_j%02d.csv", probeJ);
    writeFile(out / probeName, [&](std::ostream& f) {
        writeLineProbeCsv(f, info, mesh, U, probeJ, steps, solver.simTime());
    });
    writeFile(out / "residuals.csv", [&](std::ostream& f) {
        writeResidualsCsv(f, info, solver.residualHistory(), solver.dtHistory());
    });
    writeFile(out / "run.json", [&](std::ostream& f) {
        writeRunJson(f, info, steps, solver.simTime(), wall, status, message);
    });

    std::printf("%s  %s/%s/%s  %dx%d  steps=%d  t=%g  wall=%.3fs  status=%s\n",
                info.caseName.c_str(), info.flux.c_str(), info.reconstruction.c_str(),
                info.timeScheme.c_str(), nx, ny, steps, solver.simTime(), wall,
                status.c_str());
    if(status == "blowup")
        std::fprintf(stderr, "euler_cli: %s\n", message.c_str());

    return status == "blowup" ? 2 : 0;
}

} // anonymous namespace

int main(int argc, char** argv)
{
    Args args;
    try
    {
        args = parseArgs(argc, argv);
    }
    catch(const std::invalid_argument& e)
    {
        std::fprintf(stderr, "euler_cli: %s\nRun euler_cli --help for the options.\n", e.what());
        return 1;
    }

    if(argc == 1 || args.count("--help"))
    {
        std::fputs(USAGE, stdout);
        return argc == 1 ? 1 : 0;
    }

    try
    {
        return run(args);
    }
    catch(const std::invalid_argument& e)
    {
        std::fprintf(stderr, "euler_cli: %s\nRun euler_cli --help for the options.\n", e.what());
        return 1;
    }
    catch(const std::exception& e)
    {
        std::fprintf(stderr, "euler_cli: %s\n", e.what());
        return 1;
    }
}

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <filesystem>
#include <chrono>
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <limits>
#include <fmt/core.h>

/* Constant */
const double s2ns = 1e9, s2ms = 1e3;
const double m2um = 1e6, m2mm = 1e3;
const double EPS = std::numeric_limits<double>::epsilon();

/* TVD Runge-Kutta coefficients */
enum {nRK = 1}; // 1st-order
const double rk_weight_old[nRK] = {0.0};
const double rk_weight_new[nRK] = {1.0};

// enum {nRK = 2}; // 2nd-order
// const double rk_weight_old[nRK] = {0.0, 0.5};
// const double rk_weight_new[nRK] = {1.0, 0.5};

// enum {nRK = 3}; // 3rd-order
// const double rk_weight_old[nRK] = {0.0, 3.0/4, 1.0/3};
// const double rk_weight_new[nRK] = {1.0, 1.0/4, 2.0/3};

/* Option */
struct AppCtx
{
    double F;                // Constant-propagation speed.     Unit: m/s
    int nNode;               // Number of grid nodes.
    int nCell;               // Number of grid cells.
    double L;                // Domain length.                  Unit: m
    double h;                // Grid spacing.                   Unit: m
    double t;                // Runtime.                        Unit: s
    double tMax;             // Maximum runtime.                Unit: s
    double dt;               // Time-step.                      Unit: s
    long n;                  // Time-marching counter.
    long nMax;               // Maximum time-marching steps.
    int nOut;                // Solution output gap.
};

/* Variable */
struct Field
{
    double u;
    double phi;
};
enum {nVar = sizeof(Field) / sizeof(double)};
std::vector<std::string> VAR_NAME;

void check_range(const std::vector<Field> &v, const AppCtx &ctx)
{
    double val_min, val_max, cval;
    int min_pos, max_pos, cpos;

    // Phi @ NODE
    val_min = val_max = v[0].phi;
    min_pos = max_pos = 0;
    for (int i = 1; i < ctx.nNode; ++i)
    {
        cval = v[i].phi;
        cpos = i;
        if (cval < val_min)
        {
            val_min = cval;
            min_pos = i;
        }
        if (cval > val_max)
        {
            val_max = cval;
            max_pos = i;
        }
    }
    fmt::print("\t{:3s}: {: 14.8f} (@{:4d}) ~ {: 14.8f} (@{:4d})\n", "phi", val_min, min_pos, val_max, max_pos);
}

/* Output */
std::string runtime_str()
{
    auto tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::ostringstream ss;
    ss << std::put_time(std::localtime(&tt), "%Y%m%d-%H%M%S");
    return ss.str();
}

void write_tecplot(const std::vector<Field> &v, const AppCtx &ctx, std::ofstream &out)
{
    static const char SEP = ' ';

    out << std::setprecision(8);
    out << std::setw(18);
    out << std::right;

    out << "TITLE = \"1D Level-Set\"" << std::endl;
    out << "VARIABLES = \"X\"";
    for (int i = 0; i < nVar; i++)
        out << ", \"" << VAR_NAME[i] << "\"";
    out << std::endl;
    out << "ZONE T=\"n=" << ctx.n << "\"" << std::endl;
    out << "I=" << ctx.nCell+ctx.nNode << std::endl;
    out << "DATAPACKING=POINT" << std::endl;
    out << "SOLUTIONTIME=" << std::scientific << ctx.t << std::endl;

    for (int i = 0; i < ctx.nNode; i++)
    {
        out << std::scientific << (ctx.h * i);
        out << std::fixed;
        out << SEP << v[i].u;
        out << SEP << v[i].phi;
        out << std::endl;
    }
}

int output_solution(const std::filesystem::path &CASE_DIR, const std::vector<Field> &VAL, const AppCtx &ctx)
{
    char data_name[32];
    snprintf(data_name, sizeof(data_name), "%ld.dat", ctx.n);
    std::filesystem::path data_path = CASE_DIR/data_name;
    if (std::filesystem::exists(data_path))
    {
        std::cerr << "Target datafile \"" << data_name << "\" already exists!" << std::endl;
        return -1;
    }
    std::ofstream f_out(data_path);
    if (f_out.fail())
    {
        std::cerr << "Failed to create output datafile \"" << data_name << "\"" << std::endl;
        return -2;
    }
    write_tecplot(VAL, ctx, f_out);
    f_out.close();
    return 0;
}

/* Main program */
int main(int argc, char *argv[])
{
    AppCtx user;
    std::vector<Field> val;     // Solution @(n+1).
    std::vector<Field> val_rk;  // Solution @(m).
    std::vector<Field> val_cur; // Solution @(n).
    std::string RUN_TAG = runtime_str();
    int cnt;
    int err;

    /* Default settings */
    user.F = 1.0;
    user.nCell = 120;
    user.nNode = user.nCell + 1;
    user.L = 1.0;
    user.h = user.L / user.nCell;
    user.t = 0.0;
    user.tMax = -1.0;
    user.dt = 0.1e-9;
    user.n = 0;
    user.nMax = 0;
    user.nOut = 10;

    /* Parse options */
    cnt = 1;
    while(cnt < argc)
    {
        if (!std::strcmp(argv[cnt], "--tag"))
        {
            RUN_TAG = argv[cnt + 1];
            cnt += 2;
        }
        else if (!std::strcmp(argv[cnt], "--F"))
        {
            user.F = std::stod(argv[cnt + 1], nullptr);
            cnt += 2;
        }
        else if (!std::strcmp(argv[cnt], "--nCell"))
        {
            user.nCell = std::stoi(argv[cnt + 1], nullptr);
            user.nNode = user.nCell + 1;
            cnt += 2;
        }
        else if (!std::strcmp(argv[cnt], "--nNode"))
        {
            user.nNode = std::stoi(argv[cnt + 1], nullptr);
            user.nCell = user.nNode - 1;
            cnt += 2;
        }
        else if (!std::strcmp(argv[cnt], "--L"))
        {
            user.L = std::stod(argv[cnt + 1], nullptr);
            cnt += 2;
        }
        else if (!std::strcmp(argv[cnt], "--tMax"))
        {
            user.tMax = std::stod(argv[cnt + 1], nullptr);
            cnt += 2;
        }
        else if (!std::strcmp(argv[cnt], "--dt"))
        {
            user.dt = std::stod(argv[cnt + 1], nullptr);
            cnt += 2;
        }
        else if (!std::strcmp(argv[cnt], "--nMax"))
        {
            user.nMax = std::stol(argv[cnt + 1], nullptr);
            cnt += 2;
        }
        else if (!std::strcmp(argv[cnt], "--nOut"))
        {
            user.nOut = std::stoi(argv[cnt + 1], nullptr);
            cnt += 2;
        }
        else
        {
            std::cerr << "Option \"" << argv[cnt] << "\" is invalid!" << std::endl;
            return -1;
        }
    }

    /* Update parameters */
    user.h = user.L / user.nCell;
    user.dt = 0.5 * user.h / user.F;
    if (user.tMax < 0)
        user.tMax = user.L / user.F;
    if (user.nMax <= 0)
        user.nMax = static_cast<long>(std::ceil(user.tMax / user.dt));
    
    /* Report */
    if (!std::filesystem::exists(RUN_TAG) && !std::filesystem::create_directory(RUN_TAG))
    {
        std::cerr << "Failed to create output folder \"" << RUN_TAG << "\"" << std::endl;
        return -3;
    }
    const std::filesystem::path CASE_DIR(RUN_TAG);
    std::cout << "Output directory set to: " << CASE_DIR << std::endl;
    std::cout << "Domain length: " << user.L * m2mm << " mm" << std::endl;
    std::cout << "Number of grid cells: " << user.nCell << std::endl;
    std::cout << "Fixed time-step: " << user.dt * s2ns << " ns" << std::endl;
    std::cout << "Constant propagation speed: " << user.F  << " m/s" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;
    const std::string INFO_FILE("info.txt");
    const std::filesystem::path INFO_PATH = CASE_DIR/INFO_FILE;
    if (!std::filesystem::exists(INFO_PATH))
    {
        std::ofstream f_out(INFO_PATH);
        if (f_out.fail())
        {
            std::cerr << "Failed to create the case information file \"" << INFO_FILE << "\"" << std::endl;
            return -4;
        }
        f_out << runtime_str()                                                  << std::endl;
        f_out << "------------------------------------------------------------" << std::endl;
        f_out << "tMax = "     << user.tMax * s2ms    << " ms"       << std::endl;
        f_out << "nMax = "     << user.nMax                          << std::endl;
        f_out << "L = "        << user.L * m2mm       << " mm"       << std::endl;
        f_out << "nNode = "    << user.nNode                         << std::endl;
        f_out << "h = "        << user.h * m2um       << " um"       << std::endl;
        f_out << "dt = "       << user.dt * s2ns      << " ns"       << std::endl;
        f_out << "F = "        << user.F              << " m/s"      << std::endl;
        const double t_E = user.L / user.F * s2ms;
        f_out << "t_Entropy = " << t_E                << " ms"       << std::endl;
        f_out << "------------------------------------------------------------" << std::endl;
        f_out << "nOut = "     << user.nOut                          << std::endl;
        f_out << "------------------------------------------------------------" << std::endl;
        f_out.close();
    }

    /* Allocate storage */
    val.resize(user.nNode + user.nCell);
    val_rk.resize(user.nNode + user.nCell);
    val_cur.resize(user.nNode + user.nCell);

    /* Initialize */
    cnt = 0;
    VAR_NAME.resize(nVar);
    VAR_NAME[cnt++] = "Velocity";
    VAR_NAME[cnt++] = "SignedDistance";
    if (cnt != nVar)
    {
        std::cerr << "Unmatched variables!" << std::endl;
        return 2;
    }

    std::cout << "Apply I.C. ..." << std::endl;
    for (int idx = 0; idx < user.nNode; idx++)
    {
        val[idx].u = user.F;
        val[idx].phi = user.h * idx - user.L * 0.5;
    }
    std::copy(val.begin(), val.end(), val_rk.begin());
    std::copy(val.begin(), val.end(), val_cur.begin());

    /* Time-Marching */
    for (user.t = 0.0, user.n = 0; user.t < user.tMax && user.n < user.nMax; user.t+=user.dt, user.n++)
    {
        err = 0;
        printf("n=%ld, dt=%gms, t=%gms:\n", user.n, user.dt*s2ms, user.t*s2ms);

        /* Diagnose */
        check_range(val_cur, user);

        /* Output */
        if (user.n%user.nOut == 0)
        {
            err = output_solution(CASE_DIR, val_cur, user);
            if (err)
                break;
        }

        /* Runge-Kutta */
        for (int m = 0; m < nRK; m++)
        {
            /* Phi @NODE, Interior */
            for (int i = 1; i < user.nNode-1; i++)
            {
                double dphidx;
                if (val_rk[i].u > 0)
                    dphidx = (val_rk[i].phi - val_rk[i-1].phi) / user.h;
                else
                    dphidx = (val_rk[i+1].phi - val_rk[i].phi) / user.h;

                const double RES_LEVELSET = -val_rk[i].u * std::abs(dphidx);
                const double NEW_LEVELSET = val_rk[i].phi + user.dt * RES_LEVELSET;
                val[i].phi = rk_weight_old[m] * val_cur[i].phi + rk_weight_new[m] * NEW_LEVELSET;
            }
            /* Phi @NODE, Boundary */
            val[0].phi = val[1].phi;
            val[user.nNode-1].phi = val[user.nNode-2].phi;

            std::copy(val.begin(), val.end(), val_rk.begin());
        }

        std::copy(val.begin(), val.end(), val_cur.begin());
    }

    return err;
}

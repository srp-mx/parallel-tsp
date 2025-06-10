#pragma once

/*Copyright (C) 2025

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/


#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <dlfcn.h>
#include <math.h>

#include "defines.h"
#include "config.h"
#include "parser.h"
#include "timing.h"
#include "reporting.h"

/**
 * Macro for dynamic symbol load and error handling. Only for internal use.
 *
 * @param Err The error char* variable.
 * @param Solver The solver directly.
 * @param Name The name for the particular API attribute to load.
 */
#define LOAD_FN(Err, Solver, Name) \
    do \
    { \
        (Solver).Name = (solver_##Name *) dlsym((Solver).Handle, "solver_" #Name); \
        Err = dlerror(); \
        if (Err) \
        { \
            INSTANT_WRITE("No se pudo cargar solver_" #Name "\n"); \
            dlclose((Solver).Handle); \
            (Solver) = {}; \
            return; \
        } \
    } while(0)

typedef void command(config *Config, char *Arg);

global_variable b32 ContinueRepl = 1;

/**
 * Prints a help prompt.
 * Implements the `command` function type.
 *
 * @param Config The configuration object. Ignored.
 * @param Arg The pointer argument. Ignored.
 */
internal void
PrintHelp(config *__restrict__ Config, char *__restrict__ Arg)
{
    INSTANT_WRITE(
        "Uso:\n"
        "    > <comando> <argumento>?\n"
        "Ejemplos:\n"
        "    > problem a280\n"
        "    > solver gpu\n"
        "    > help\n"
        "Comandos:\n"
        "    help: Muestra este diálogo\n"
        "    problem <nombre de archivo>: Carga un archivo de problema\n"
        "    solver <nombre del solucionador>: Asigna el algoritmo del solucionador\n"
        "    iterations <entero positivo>: Asigna el número de iteraciones por ejecución\n"
        "    executions <entero positivo>: Asigna el número de ejecuciones del algoritmo\n"
        "    cutoff <número positivo>: Asigna el valor suficiente de costo de la solución\n"
        "    parallelism <entero positivo o -1>: Asigna el nivel de paralelismo, sujeto al solucionador\n"
        "    config: Imprime los datos de la configuración actual\n"
        "    exit: Cierra el programa\n"
        "    run: Ejecuta el programa con la configuración dada\n"
    );
}

/**
 * Loads a problem file into the configuration.
 * Implements the `command` function type.
 *
 * @param Config The configuration object.
 * @param Arg The pointer argument.
 */
internal void
LoadProblem(config *__restrict__ Config, char *__restrict__ Arg)
{
    // Wrap the argument around the directory and extension
    char FileName[MAX_FILENAME];
    FileName[0] = '.';
    FileName[1] = '/';
    FileName[2] = 't';
    FileName[3] = 's';
    FileName[4] = 'p';
    FileName[5] = '/';
    size_t ArgLen = strnlen(Arg, MAX_FILENAME-30);
    memcpy(FileName+6, Arg, ArgLen);
    FileName[6+ArgLen] = '.';
    FileName[6+ArgLen+1] = 't';
    FileName[6+ArgLen+2] = 's';
    FileName[6+ArgLen+3] = 'p';
    FileName[6+ArgLen+4] = 0;

    // Search and load the file
    struct stat Stat;
    if (Config->ProblemFileSize > 0)
    {
        munmap(Config->ProblemSource, Config->ProblemFileSize);
        Config->ProblemFileSize = 0;
        Config->ProblemSource = 0;
        Config->Name = 0;
        Config->Comment = 0;
    }
    if (Config->Tsp)
    {
        munmap(Config->Tsp, sizeof(tsp_instance) + sizeof(v2)*Config->Tsp->N);
        Config->Tsp = 0;
    }
    if (stat(FileName, &Stat))
    {
        INSTANT_WRITE("No existe el archivo.\n");
        return;
    }
    Config->ProblemFileSize = Stat.st_size;
    int FD = open(FileName, 0);
    Config->ProblemSource = (char *)mmap(0, Config->ProblemFileSize,
            PROT_READ | PROT_WRITE,
            MAP_PRIVATE,
            FD, 0);
    close(FD);

    // Parse the file
    parsed_tsp Parsed = ParseTsp(Config->ProblemSource);
    if (!Parsed.Ok)
    {
        INSTANT_WRITE("El contenido del archivo no es válido.\n");
        if (Parsed.Tsp)
        {
            munmap(Parsed.Tsp, sizeof(tsp_instance) + sizeof(v2)*Parsed.N);
        }
        return;
    }
    Config->Tsp = Parsed.Tsp;
    Config->Name = Parsed.Name;
    Config->Comment = Parsed.Comment;
}

/**
 * Loads a solver dynamic library file into the program.
 * Implements the `command` function type.
 *
 * @param Config The configuration object.
 * @param Arg The pointer argument.
 */
internal void
LoadSolver(config *__restrict__ Config, char *__restrict__ Arg)
{
    if (Config->Solver.Handle)
    {
        Config->Solver.Unload();
        dlclose(Config->Solver.Handle);
        Config->Solver = {};
    }

    char LibPath[MAX_FILENAME];
    LibPath[0] = '.';
    LibPath[1] = '/';
    size_t ArgLen = strnlen(Arg, MAX_FILENAME-30);
    memcpy(LibPath+2, Arg, ArgLen);
    LibPath[2+ArgLen] = '.';
    LibPath[2+ArgLen+1] = 's';
    LibPath[2+ArgLen+2] = 'o';
    LibPath[2+ArgLen+3] = 0;

    // Open the dynamic library
    Config->Solver.Handle = dlopen(LibPath, RTLD_LOCAL | RTLD_NOW);
    if (!Config->Solver.Handle)
    {
        INSTANT_WRITE("No se pudo abrir el solucionador:\n>>>>>>>\n");
        const char *Err = dlerror();
        write(1, Err, strlen(Err));
        INSTANT_WRITE("\n<<<<<<<\n");
        return;
    }

    // Clear existing errors
    dlerror();

    // Load the compatibility version function pointer and verify compatibility
    const char *Err = 0;
    LOAD_FN(Err, Config->Solver, Compatibility);
    if (Config->Solver.Compatibility() != MAIN_MAJOR_VERSION)
    {
        dlclose(Config->Solver.Handle);
        Config->Solver = {};
        INSTANT_WRITE("El solucionador no es compatible con esta versión del programa.\n");
        return;
    }

    // Load all remaining function pointers
    LOAD_FN(Err, Config->Solver, Version);
    LOAD_FN(Err, Config->Solver, Name);
    LOAD_FN(Err, Config->Solver, Description);
    LOAD_FN(Err, Config->Solver, Setup);
    LOAD_FN(Err, Config->Solver, Unload);
    LOAD_FN(Err, Config->Solver, Solve);

    // Once everything is in order, run the library setup
    if (!Config->Solver.Setup())
    {
        dlclose(Config->Solver.Handle);
        Config->Solver = {};
        INSTANT_WRITE("Hubo un error inicializando el solucionador.\n");
        return;
    }
}

/**
 * Sets the number of iterations of each execution of the program.
 * Implements the `command` function type.
 *
 * @param Config The configuration object.
 * @param Arg The pointer argument to the number to store.
 */
internal void
SetIterations(config *__restrict__ Config, char *__restrict__ Arg)
{
    char *End;
    Config->Iterations = strtoull(Arg, &End, 10);
    if (!Config->Iterations || Config->Iterations >> 63)
    {
        INSTANT_WRITE("Se esperaba un entero positivo, lo cual no se leyó.\n");
        Config->Iterations = 0;
    }
}

/**
 * Sets the number of executions of the solver with this configuration.
 * Implements the `command` function type.
 *
 * @param Config The configuration object.
 * @param Arg The pointer argument to the number to store.
 */
internal void
SetExecutions(config *__restrict__ Config, char *__restrict__ Arg)
{
    char *End;
    Config->Executions = strtoull(Arg, &End, 10);
    if (!Config->Executions || Config->Executions >> 63)
    {
        INSTANT_WRITE("Se esperaba un entero positivo, lo cual no se leyó.\n");
        Config->Executions = 0;
    }
}

/**
 * Sets the upper bound on the target cost. If the solver's cost reaches or
 * improves upon this, it should exit and return the solution.
 * Implements the `command` function type.
 *
 * @param Config The configuration object.
 * @param Arg The pointer argument to the value to store.
 */
internal void
SetCutoff(config *__restrict__ Config, char *__restrict__ Arg)
{
    char *End;
    Config->Cutoff = strtof(Arg, &End);
    if (Config->Cutoff <= 0.0f)
    {
        INSTANT_WRITE("Se esperaba un número positivo, lo cual no se leyó.\n");
        Config->Cutoff = 0.0f;
    }
}

/**
 * Sets the level of parallelism to use, whose meaning depends on the solver.
 *
 * @param Config The configuration object.
 * @param Arg The pointer argument to the value to store.
 */
internal void
SetParallelism(config *__restrict__ Config, char *__restrict__ Arg)
{
    char *End;
    Config->Parallelism = (i32)strtol(Arg, &End, 10);
    if (!Config->Parallelism || Config->Parallelism < -1)
    {
        INSTANT_WRITE("Se esperaba un entero positivo o -1, lo cual no se leyó.\n");
        Config->Iterations = 0;
    }

    if (Config->Parallelism == -1) Config->Parallelism = 0;
}

/**
 * Prints the state of the execution configuration.
 * Implements the `command` function type.
 *
 * @param Config The configuration object.
 * @param Arg The pointer argument. Ignored.
 */
internal void
PrintConfig(config *__restrict__ Config, char *__restrict__ Arg)
{
    char *ProblemSource = Config->ProblemSource;
    u64 Iterations = Config->Iterations;
    u64 Executions = Config->Executions;
    solver *Solver = &(Config->Solver);
    if (!ProblemSource && Iterations <= 1 && Executions <= 1 && !Solver->Handle)
    {
        INSTANT_WRITE("No se ha cargado ninguna información de ejecución.\n");
        return;
    }

    size_t Bytes = ConfigToStr(Config, 0);
    char OutBuff[Bytes] = {};
    Bytes = ConfigToStr(Config, OutBuff);
    write(1, OutBuff, Bytes);
}

/**
 * Indicates the REPL must exit.
 *
 * @param Config Ignored.
 * @param Arg Ignored.
 */
internal void
Exit(config *__restrict__ Config, char *__restrict__ Arg)
{
    ContinueRepl = 0;
    INSTANT_WRITE("Adiós :)\n");
}

/**
 * Executes the program according to the configuration and logs results.
 *
 * @param Config Pointer to the configuration.
 * @param Arg Ignored.
 */
internal void
Run(config *__restrict__ Config, char *__restrict__ Arg)
{
    if (!Config->Tsp)
    {
        INSTANT_WRITE("No se ha cargado una instancia del problema. Usa help.\n");
        return;
    }
    if (!Config->Solver.Handle)
    {
        INSTANT_WRITE("No se ha cargado un solucionador. Usa help.\n");
        return;
    }
    if (Config->Iterations < 1)
    {
        INSTANT_WRITE("El número de iteraciones debe ser positivo.\n");
        return;
    }
    if (Config->Executions < 1)
    {
        INSTANT_WRITE("El número de ejecuciones debe ser positivo.\n");
        return;
    }

    size_t ProblemBytes = sizeof(tsp_instance) + sizeof(v2)*Config->Tsp->N;
    tsp_instance *ProblemCopy = (tsp_instance*)alloca(ProblemBytes);
    r32 SolutionCosts[Config->Executions] = {};
    r32 BestValue = INFINITY;
    i32 BestRoute[Config->Tsp->N] = {};
    i32 BestIdx = 0;
    r64 Timings[Config->Executions] = {};
    u64 TimingsRdtsc[Config->Executions] = {};
    u64 Iterations[Config->Executions] = {};
    r32 Cutoff = Config->Cutoff;
    r32 HitPercent = 0.0f;
    i32 Parallelism = Config->Parallelism;

    for (u64 I = 0; I < Config->Executions; I++)
    {
        Iterations[I] = Config->Iterations;
    }

    for (u64 Execution = 0; Execution < Config->Executions; Execution++)
    {
        memcpy(ProblemCopy, Config->Tsp, ProblemBytes);
        i32 Route[Config->Tsp->N] = {};

        b32 Ok;
        MEASURE_TIME(Time,
            Ok = Config->Solver.Solve(ProblemCopy, Route, Iterations+Execution,
                    Cutoff, Parallelism);
        );
        Timings[Execution] = fsecs_Time;
        TimingsRdtsc[Execution] = rdtsc_Time;

        if (!Ok)
        {
            INSTANT_WRITE("El solucionador reportó un problema.\n");
            return;
        }

        // Verify all vertices were explored exactly once
        u8 Explored[Config->Tsp->N] = {};
        for (i32 I = 0; I < Config->Tsp->N; I++)
        {
            if (Route[I] < 0 || Route[I] >= Config->Tsp->N)
            {
                INSTANT_WRITE("El solucionador dio un índice fuera de rango.\n");
                return;
            }
            if (Explored[Route[I]])
            {
                INSTANT_WRITE("El solucionador dio una solución que repite "
                        "vértices.\n");
                for (i32 i = 0; i < Config->Tsp->N; i++)
                {
                    printf("%d ", Route[i]);
                }
                printf("\n");
                return;
            }
            Explored[Route[I]] = 1;
        }
        for (i32 I = 0; I < Config->Tsp->N; I++)
        {
            if (!Explored[Route[I]])
            {
                INSTANT_WRITE("El solucionador dio una solución que no pasa "
                        "por todos los vértices.\n");
                return;
            }
        }

        // Obtain the cost of the route
        r32 RouteCost = 0.0f;
        for (i32 I = 0; I < Config->Tsp->N; I++)
        {
            v2 U = Config->Tsp->Coords[Route[I]];
            v2 V = Config->Tsp->Coords[Route[(I+1)%Config->Tsp->N]];
            RouteCost += Dist(U, V);
        }

        // Store stats
        SolutionCosts[Execution] = RouteCost;
        if (RouteCost <= Cutoff)
        {
            HitPercent += 1.0f;
        }

        // Set best if this is it
        if (RouteCost < BestValue)
        {
            BestIdx = Execution;
            BestValue = RouteCost;
            memcpy(BestRoute, Route, sizeof(i32)*Config->Tsp->N);
        }
    }
    HitPercent /= Config->Executions;

    report_data ReportData = {};
    ReportData.Problem = Config->Name;
    ReportData.Solver = Config->Solver.Name();
    ReportData.Costs = SolutionCosts;
    ReportData.BestRoute = BestRoute;
    ReportData.Times = Timings;
    ReportData.Cycles = TimingsRdtsc;
    ReportData.Iterations = Iterations;
    ReportData.BestCost = BestValue;
    ReportData.BestCostIdx = BestIdx;
    ReportData.N = Config->Tsp->N;
    ReportData.Execs = Config->Executions;
    ReportData.HitPercent = HitPercent;
    ReportData.Cutoff = Cutoff;
    ReportData.Parallelism = Config->Parallelism;

    // Create the report
    char Report[ReportToStr(&ReportData, 0)];
    size_t ReportLen = ReportToStr(&ReportData, Report);

    // Write it into console
    write(1, Report, ReportLen);

    // Name the outgoing file datetime%problem%solver.json
    char ReportFilePath[DateTimeStr(0) + strlen(ReportData.Problem)
        + strlen(ReportData.Solver) + sizeof("./logs/.json") + 1] = {};
    char Datetime[DateTimeStr(0)] = {}; DateTimeStr(Datetime);
    sprintf(ReportFilePath, "./logs/%s%%%s%%%s.json",
            Datetime, ReportData.Problem, ReportData.Solver);

    // Cleanup filename to remove spaces
    for (char *C = ReportFilePath; *C; C++)
    {
        if (*C == ' ') *C = '_';
    }

    // Save the report into a new file
    i32 Fd = open(ReportFilePath,
            O_CREAT | O_WRONLY | O_TRUNC,
            0644);
    write(Fd, Report, ReportLen);
    close(Fd);
}

// Number of commands
#define CMD_COUNT 10

// List of command names
global_variable const char *COMMAND_NAMES[CMD_COUNT] = {
    "help",
    "problem",
    "solver",
    "iterations",
    "executions",
    "cutoff",
    "parallelism",
    "config",
    "exit",
    "run"
};

// List of command name lengths
global_variable const size_t COMMAND_NAME_LENS[CMD_COUNT] = {
    size_t(strlen("help")),
    size_t(strlen("problem")),
    size_t(strlen("solver")),
    size_t(strlen("iterations")),
    size_t(strlen("executons")),
    size_t(strlen("cutoff")),
    size_t(strlen("parallelism")),
    size_t(strlen("config")),
    size_t(strlen("exit")),
    size_t(strlen("run"))
};

// List of command functions, lined up with `COMMAND_NAMES`
global_variable const command *COMMANDS[CMD_COUNT] = {
    PrintHelp,
    LoadProblem,
    LoadSolver,
    SetIterations,
    SetExecutions,
    SetCutoff,
    SetParallelism,
    PrintConfig,
    Exit,
    Run
};

#undef LOAD_FN

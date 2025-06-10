#include "solver.h"

#include "solve.cuh"

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

#define ERR_DEVCNT "Couldn't check CUDA compatibility.\n"
#define ERR_NOCUDA "This machine has no CUDA capable devices.\n"
#define WARN_EXTRAGPUS "Warning, this machine has multiple CUDA devices.\n" \
                       "    Only device 0 will be used.\n"

global_variable cudaDeviceProp DeviceProperties;

/**
 * Returns the version of the solver.
 */
u64
solver_Version()
{
    return 2;
}

/**
 * Returns the version of the main program the solver is compatible with.
 */
u64
solver_Compatibility()
{
    return 2;
}

/**
 * Returns the name of the solver.
 */
const char *
solver_Name()
{
    return "GPU Parallel Genetic Solver";
}

/**
 * Returns the description of the solver.
 */
const char *
solver_Description()
{
    return "Closely follows the elite island model from the paper at "
        "http://difu100cia.uaz.edu.mx/index.php/difuciencia/article/view/145";
}

/**
 * Checks if there is a CUDA device and loads its data.
 *
 * @return 1 (everything ok) or 0 (error).
 */
b32
solver_Setup()
{
    i32 DeviceCount = 0;
    cudaError_t Err = cudaGetDeviceCount(&DeviceCount);

    if (Err != cudaSuccess)
    {
        IGNORE_RESULT(write(1, ERR_DEVCNT, sizeof(ERR_DEVCNT)));
        return 0;
    }

    if (DeviceCount == 0)
    {
        IGNORE_RESULT(write(1, ERR_NOCUDA, sizeof(ERR_NOCUDA)));
        return 0;
    }
    else if (DeviceCount > 1)
    {
        IGNORE_RESULT(write(1, WARN_EXTRAGPUS, sizeof(WARN_EXTRAGPUS)));
        return 0;
    }

    cudaSetDevice(0);
    cudaGetDeviceProperties(&DeviceProperties, 0);

    return 1;
}

/**
 * Does any necessary library unloading procedures.
 */
void
solver_Unload() {}

/**
 * Executes the euclidean TSP solver.
 *
 * @param Tsp A pointer to the TSP instance to be read from.
 * @param out_Permutation An array to be filled in with the best solution.
 * @param Iterations Has the maximum number of iterations and the ones executed
 *                   are written back into it.
 * @param Cutoff If this value is reached or passed, we stop immediately.
 * @param Parallelism Number of islands (thread blocks), if zero lets the
 *                    solver decide.
 * 
 * @return 1 if everything is ok, or 0 otherwise.
 */
b32
solver_Solve(tsp_instance *__restrict__ Tsp,
             i32 *__restrict__ out_Permutation,
             u64 *__restrict__ Iterations,
             r32 Cutoff,
             i32 Parallelism)
{
    return Main(Tsp, out_Permutation, Iterations, &DeviceProperties, Cutoff, Parallelism);
}

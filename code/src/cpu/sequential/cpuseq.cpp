#include "solver.h"

#include "solve.cpp"

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

/**
 * Returns the version of the solver.
 */
u64
solver_Version()
{
    return 1;
}

/**
 * Returns the version of the main program the solver is compatible with.
 */
u64
solver_Compatibility()
{
    return 1;
}

/**
 * Returns the name of the solver.
 */
const char *
solver_Name()
{
    return "CPU Sequential Genetic Solver";
}

/**
 * Returns the description of the solver.
 */
const char *
solver_Description()
{
    return "Follows the concept explained in the dissertation at\n"
        "http://132.248.9.195/ptd2023/septiembre/0846778/Index.html";
}

/**
 * Checks if there is a CUDA device and loads its data.
 *
 * @return 1 (everything ok) or 0 (error).
 */
b32
solver_Setup()
{
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
 * @param out_Permutation An array with N spaces, to be 0 through N-1.
 * @param Iterations A single iteration will be executed, so it will write 1.
 * 
 * @return 1 (everything ok).
 */
b32
solver_Solve(tsp_instance *__restrict__ Tsp,
             i32 *__restrict__ out_Permutation,
             u64 *__restrict__ Iterations,
             r32 Cutoff)
{
    return Main(Tsp, out_Permutation, Iterations, Cutoff);
}

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


// The following are properly documented in src/include/solver.h
typedef u64 solver_Version();
typedef u64 solver_Compatibility();
typedef const char *solver_Name();
typedef const char *solver_Description();
typedef b32 solver_Setup();
typedef void solver_Unload();
typedef b32 solver_Solve(tsp_instance *Tsp,
                         i32 *out_Permutation,
                         u64 *Iterations,
                         r32 Cutoff,
                         i32 Parallelism);

/**
 * Structure to hold all solver function pointers.
 */
struct solver
{
    void *Handle;
    solver_Version *Version;
    solver_Compatibility *Compatibility;
    solver_Name *Name;
    solver_Description *Description;
    solver_Setup *Setup;
    solver_Unload *Unload;
    solver_Solve *Solve;
};

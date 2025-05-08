#pragma once

#include "defines.h"
#include "tsp.h"

#define EXPORT __attribute__((visibility("default"))) 

/**
 * Returns the version of the solver.
 */
extern "C" EXPORT u64
solver_Version();

/**
 * Returns the version of the main program the solver is compatible with.
 */
extern "C" EXPORT u64
solver_Compatibility();

/**
 * Returns the name of the solver.
 */
extern "C" EXPORT const char *
solver_Name();

/**
 * Returns the description of the solver.
 */
extern "C" EXPORT const char *
solver_Description();

/**
 * Does any necessary library setup and indicates if it was successful.
 *
 * @return 1 if successful.
 */
extern "C" EXPORT b32
solver_Setup();

/**
 * Does any necessary library unloading procedures.
 */
extern "C" EXPORT void
solver_Unload();

/**
 * Gives an approximate solution to the euclidean TSP.
 *
 * @param Tsp A pointer to the TSP instance to be read from, and optionally
 *            written into.
 * @param out_Permutation The node permutation given as a solution. Any changes
 *                        to the Tsp variable will be ignored. It is already
 *                        allocated for you and has as many entries as there
 *                        are nodes in the problem instance definition.
 * @param Iterations Pointer to the number of iterations to run the solver,
 *                   which should be written back to with the number of
 *                   iterations actually executed.
 * @param Cutoff The minimum cost deemed as a good solution. If a solution with
 *               this value or less is found, the solver HAS to terminate and
 *               return it.
 * 
 * @return 1 if everything went ok.
 */
extern "C" EXPORT b32
solver_Solve(tsp_instance *__restrict__ Tsp,
             i32 *__restrict__ out_Permutation,
             u64 *__restrict__ Iterations,
             r32 Cutoff);

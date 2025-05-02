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

#include "tsp.h"
#include "api.h"

#define MAX_FILENAME 512
#define NO_PROBLEM_NAME "<Sin Nombre>"
#define NO_PROBLEM_COMMENT "<Sin Comentarios>"
#define FMTSTR_SOLVER "SOLUCIONADOR: %s [VERSIÓN: %lu]\n    %s\n"
#define FMTSTR_PROBLEM "PROBLEMA: %s\n    %s\n"
#define FMTSTR_ITERS "ITERACIONES: %lu\n"
#define FMTSTR_EXECS "EJECUCIONES: %lu\n"

/**
 * Contains the execution configuration.
 */
struct config
{
    size_t ProblemFileSize;
    char *ProblemSource;
    tsp_instance *Tsp;
    char *Name;
    char *Comment;
    solver Solver;
    u64 Iterations, Executions;
};

/**
 * Writes the string representation of the current configuration.
 * If no buffer is passed in (null pointer), then nothing is written.
 *
 * @param Config The configuration object to convert.
 * @param Buffer Where to write the string, ignored if null.
 *
 * @return How many bytes were used (excluding the terminator), or if Buffer is
 *         null, how many should be allocated for the string to fit in the
 *         worst case.
 */
internal size_t
ConfigToStr(config *Config, char *Buffer)
{
    size_t ProblemNameLen =
        Config->Name ? strlen(Config->Name) : sizeof(NO_PROBLEM_NAME);
    size_t ProblemCommentLen =
        Config->Comment ? strlen(Config->Comment) : sizeof(NO_PROBLEM_COMMENT);
    const char *SolverName = "";
    const char *SolverDesc = "";
    u64 SolverVers = 0;
    if (Config->Solver.Handle)
    {
        SolverName = Config->Solver.Name();
        SolverDesc = Config->Solver.Description();
        SolverVers = Config->Solver.Version();
    }

    // If there is no buffer, just estimate an upper bound on the bytes used.
    if (!Buffer)
    {
        return sizeof(FMTSTR_SOLVER) + sizeof(FMTSTR_PROBLEM)
            + sizeof(FMTSTR_ITERS) + sizeof(FMTSTR_EXECS) + ProblemNameLen
            + ProblemCommentLen + 20*3 + strlen(SolverName)
            + strlen(SolverDesc) + 1;
    }

    // If there is a buffer, just sprintf it.
    size_t Written = 0;
    if (Config->Solver.Handle)
    {
        Written += sprintf(Buffer+Written, FMTSTR_SOLVER,
                SolverName, SolverVers, SolverDesc);
    }
    if (Config->Tsp)
    {
        const char *Name = Config->Name ? Config->Name : NO_PROBLEM_NAME;
        const char *Comm = Config->Comment ? Config->Comment : NO_PROBLEM_COMMENT;
        Written += sprintf(Buffer+Written, FMTSTR_PROBLEM, Name, Comm);
    }
    if (Config->Executions != 1)
    {
        Written += sprintf(Buffer+Written, FMTSTR_EXECS, Config->Executions);
    }
    if (Config->Iterations != 1)
    {
        Written += sprintf(Buffer+Written, FMTSTR_ITERS, Config->Iterations);
    }

    return Written;
}

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

#include "timing.h"
#include "utils.h"

/**
 * NOTE(srp): JSON specification for the report
 *
 * {
 *     "metadata": {
 *         "datetime": "...",
 *         "problem": "...",
 *         "solver": "...",
 *         "n": 123,
 *         "executions": 123
 *     },
 *     "best_solution": {
 *         "route": [ 1, 2, 3, ... ],
 *         "cost": 123.456,
 *         "idx": 7
 *     },
 *     "stats": {
 *         "cost": {
 *             "avg": 12.34,
 *             "stdev": 12.34,
 *             "min": 12.34,
 *             "max": 12.34
 *         },
 *         "seconds": {
 *             "avg": 12.34,
 *             "stdev": 12.34,
 *             "min": 12.34,
 *             "max": 12.34
 *         },
 *         "cycles": {
 *             "avg": 12.34,
 *             "stdev": 12.34,
 *             "min": 1234,
 *             "max": 1234
 *         },
 *         "iterations": {
 *             "avg": 12.34,
 *             "stdev": 12.34,
 *             "min": 1234,
 *             "max": 1234
 *         }
 *     },
 *     "raw_metrics": {
 *         "costs": [ 1.23, 1.23, ... ],
 *         "seconds": [ 1.23, 1.23, ... ],
 *         "cycles": [ 123, 123, ... ],
 *         "iterations": [ 1.23, 1.23, ... ]
 *     }
 * }
 */

#define FMTSTR_METADATA \
    "{\n" \
    "    \"metadata\": {\n" \
    "        \"datetime\": \"%s\",\n" \
    "        \"problem\": \"%s\",\n" \
    "        \"solver\": \"%s\",\n" \
    "        \"n\": %d,\n" \
    "        \"executions\": %d\n" \
    "    },\n"

#define FMTSTR_BESTSOL \
    "    \"best_solution\": {\n" \
    "        \"route\": [ %s ],\n" \
    "        \"cost\": %e,\n" \
    "        \"idx\": %d\n" \
    "    },\n"

#define FMTSTR_STATSGENERIC(FmtElement) \
    ": {\n" \
    "            \"avg\": %e,\n" \
    "            \"stdev\": %e,\n" \
    "            \"min\": " FmtElement ",\n" \
    "            \"max\": " FmtElement "\n" \
    "        }"

#define FMTSTR_STATS \
    "    \"stats\": {\n" \
    "        \"cost\"" FMTSTR_STATSGENERIC("%e") ",\n" \
    "        \"seconds\"" FMTSTR_STATSGENERIC("%e") ",\n" \
    "        \"cycles\"" FMTSTR_STATSGENERIC("%lu") ",\n" \
    "        \"iterations\"" FMTSTR_STATSGENERIC("%lu") "\n" \
    "    },\n"

#define FMTSTR_RAWMETRICS \
    "    \"raw_metrics\": {\n" \
    "        \"costs\": [ %s ],\n" \
    "        \"seconds\": [ %s ],\n" \
    "        \"cycles\": [ %s ],\n" \
    "        \"iterations\": [ %s ]\n" \
    "    }\n" \
    "}\n"

#define FMTSTR_REPORT \
    ( FMTSTR_METADATA FMTSTR_BESTSOL FMTSTR_STATS FMTSTR_RAWMETRICS )

/**
 * Structure with information to include in the report.
 */
struct report_data
{
    char *Problem;
    const char *Solver;
    r32 *Costs;
    i32 *BestRoute;
    r64 *Times;
    u64 *Cycles;
    u64 *Iterations;
    r32 BestCost;
    i32 BestCostIdx;
    i32 N;
    i32 Execs;
};

/**
 * Writes a report into the provided string buffer and returns how many bytes
 * were used. If Buffer is null, it only returns an upper bound estimate of
 * how many bytes are required to store the report as a string.
 *
 * @param Data The report data, completely filled in.
 * @param Buffer Null or the address to write the string into.
 *
 * @return The amount of bytes used, or the amount required if Buffer is null.
 */
internal size_t
ReportToStr(report_data *Data, char *Buffer)
{
    if (!Buffer)
    {
        size_t Written = 0;
        Written += sizeof(FMTSTR_REPORT)+1;
        Written += DateTimeStr(0) + strlen(Data->Problem)
            + strlen(Data->Solver) + 20*2; // metadata values
        Written += 20*Data->N + 25 + 20; // best_solution values
        Written += 20*2*4 + 20*2*2 + 25*2*2; // stats values
        Written += 2*25*Data->Execs + 2*20*Data->Execs; // raw_metrics values
        return Written;
    }

    // First, obtain everything that will be formatted as a string.
    char Datetime[DateTimeStr(0)];
    DateTimeStr(Datetime);
    char RouteStr[20*Data->N];
    char *S0 = RouteStr;
    for (i32 I = 0; I < Data->N; I++)
    {
        if (I < Data->N - 1)
        {
            S0 += sprintf(S0, "%d, ", (Data->BestRoute[I]) + 1);
        }
        else
        {
            S0 += sprintf(S0, "%d", (Data->BestRoute[I]) + 1);
        }
    }
    char CostsStr[25*Data->Execs]; char *S1 = CostsStr;
    char SecondsStr[25*Data->Execs]; char *S2 = SecondsStr;
    char CyclesStr[20*Data->Execs]; char *S3 = CyclesStr;
    char IterationsStr[20*Data->Execs]; char *S4 = IterationsStr;
    for (i32 I = 0; I < Data->Execs; I++)
    {
        if (I < Data->Execs - 1)
        {
            S1 += sprintf(S1, "%e, ", Data->Costs[I]);
            S2 += sprintf(S2, "%e, ", Data->Times[I]);
            S3 += sprintf(S3, "%lu, ", Data->Cycles[I]);
            S4 += sprintf(S4, "%lu, ", Data->Iterations[I]);
        }
        else
        {
            S1 += sprintf(S1, "%e", Data->Costs[I]);
            S2 += sprintf(S2, "%e", Data->Times[I]);
            S3 += sprintf(S3, "%lu", Data->Cycles[I]);
            S4 += sprintf(S4, "%lu", Data->Iterations[I]);
        }
    }

    // Do the string formatting
    i32 E = Data->Execs;
    return sprintf(Buffer, FMTSTR_REPORT,
        // "metadata"
        Datetime, Data->Problem, Data->Solver, Data->N, Data->Execs,
        // "best_solution"
        RouteStr, Data->BestCost, Data->BestCostIdx,
        // "stats.cost"
        r32_Mean(E, Data->Costs), r32_Stdev(E, Data->Costs),
        r32_Min(E, Data->Costs), r32_Max(E, Data->Costs),
        // "stats.seconds"
        r64_Mean(E, Data->Times), r64_Stdev(E, Data->Times),
        r64_Min(E, Data->Times), r64_Max(E, Data->Times),
        // "stats.cycles"
        u64_Mean(E, Data->Cycles), u64_Stdev(E, Data->Cycles),
        u64_Min(E, Data->Cycles), u64_Max(E, Data->Cycles),
        // "stats.iterations"
        u64_Mean(E, Data->Iterations), u64_Stdev(E, Data->Iterations),
        u64_Min(E, Data->Iterations), u64_Max(E, Data->Iterations),
        // "raw_metrics"
        CostsStr, SecondsStr, CyclesStr, IterationsStr);
}

#undef FMTSTR_METADATA
#undef FMTSTR_BESTSOL
#undef FMTSTR_STATSGENERIC
#undef FMTSTR_STATS
#undef FMTSTR_RAWMETRICS
#undef FMTSTR_REPORT

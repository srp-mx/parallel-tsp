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

#include <time.h>
#include <immintrin.h>

/**
 * Macro to wrap around an expression and measure its timing.
 *
 * @param MeasurementName Name of the resulting variables with the prefix
 *                        rdtsc_ and fsecs_ for the time as CPU cycles and
 *                        seconds (double precision) respectively.
 * @param Expression The expression to measure.
 */
#define MEASURE_TIME(MeasurementName, Expression) \
    u64 rdtsc_##MeasurementName; \
    r64 fsecs_##MeasurementName; \
    do \
    { \
        u64 IniRdtsc_, EndRdtsc_; \
        timespec IniTspec_, EndTspec_; \
        clock_gettime(CLOCK_MONOTONIC, &IniTspec_); \
        IniRdtsc_ = __rdtsc(); \
        do { Expression; } while(0); \
        EndRdtsc_ = __rdtsc(); \
        clock_gettime(CLOCK_MONOTONIC, &EndTspec_); \
        rdtsc_##MeasurementName = EndRdtsc_ - IniRdtsc_; \
        fsecs_##MeasurementName = \
            (EndTspec_.tv_sec + ((r64)EndTspec_.tv_nsec)/1e9) - \
            (IniTspec_.tv_sec + ((r64)IniTspec_.tv_nsec)/1e9); \
    } while(0);

/**
 * Writes the datetime into a buffer and returns how much bytes were used.
 *
 * @param Buffer The buffer to write into.
 *
 * @return The number of bytes used/to use.
 */
internal size_t
DateTimeStr(char *Buffer)
{
    constexpr size_t FMT_LEN = sizeof("yyyy-mm-ddThh_mm_ss.sssZ");
    if (!Buffer)
    {
        return FMT_LEN+1;
    }

    timespec Ts;
    tm Tm;
    clock_gettime(CLOCK_REALTIME, &Ts);
    localtime_r(&Ts.tv_sec, &Tm);
    return snprintf(Buffer, FMT_LEN + 1,
            "%04d-%02d-%02dT%02d_%02d_%02d.%03ldZ",
            Tm.tm_year + 1900, Tm.tm_mon + 1, Tm.tm_mday,
            Tm.tm_hour, Tm.tm_min, Tm.tm_sec,
            Ts.tv_nsec / 1000000);
}

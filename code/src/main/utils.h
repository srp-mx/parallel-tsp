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

#include <math.h>

/**
 * Gets the maximum between two size_t instances.
 *
 * @param A The first element.
 * @param B The second element.
 *
 * @return The largest between A and B
 */
inline internal size_t
Max(size_t A, size_t B)
{
    return A <= B ? B : A;
}

/**
 * Macro to define mean functions with numeric types.
 */
#define MEAN_FN(Type) \
    internal r64 \
    Type##_Mean(i32 N, Type *Data) \
    { \
        r64 Result = 0; \
        for (i32 I = 0; I < N; I++) \
        { \
            Result += (r64)Data[I]; \
        } \
        return Result/(r64)N; \
    }

/**
 * Macro to define sample standard deviation functions with numeric types.
 */
#define STDEV_FN(Type) \
    internal r64 \
    Type##_Stdev(i32 N, Type *Data) \
    { \
        r64 Mu = Type##_Mean(N, Data); \
        r64 Result = 0; \
        for (i32 I = 0; I < N; I++) \
        { \
            r64 Xi = (r64)Data[I] - Mu; \
            Result += Xi*Xi; \
        } \
        return sqrt(Result/(r64)(N-1)); \
    }

/**
 * Macro to define min functions with a particular type.
 */
#define MIN_FN(Type) \
    internal Type \
    Type##_Min(i32 N, Type *Data) \
    { \
        Type Min = Data[0]; \
        for (i32 I = 1; I < N; I++) \
        { \
            if (Data[I] < Min) Min = Data[I]; \
        } \
        return Min; \
    }

/**
 * Macro to define max functions with a particular type.
 */
#define MAX_FN(Type) \
    internal Type \
    Type##_Max(i32 N, Type *Data) \
    { \
        Type Max = Data[0]; \
        for (i32 I = 1; I < N; I++) \
        { \
            if (Data[I] > Max) Max = Data[I]; \
        } \
        return Max; \
    }

// Definition for multiple statistics function with different numeric types.
MEAN_FN(u64)
MEAN_FN(r64)
MEAN_FN(r32)
STDEV_FN(u64)
STDEV_FN(r64)
STDEV_FN(r32)
MIN_FN(u64)
MIN_FN(r64)
MIN_FN(r32)
MAX_FN(u64)
MAX_FN(r64)
MAX_FN(r32)

#undef MEAN_FN
#undef STDEV_FN
#undef MIN_FN
#undef MAX_FN

/**
 * Returns the euclidean distance between two points (L2 norm).
 * 
 * @param U The first point.
 * @param V The second point.
 *
 * @return The euclidean distance between U and V.
 */
inline internal r32
Dist(v2 U, v2 V)
{
    r32 Dx = U.X - V.X;
    r32 Dy = U.Y - V.Y;
    return sqrtf(Dx*Dx + Dy*Dy);
}

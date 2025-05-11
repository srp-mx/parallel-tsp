#pragma once

#include "defines.h"

#include <unistd.h>
#include <stdio.h>
#include <assert.h>
#include <immintrin.h>

#include <curand_kernel.h>

#define STRFMT_GLOBUSE "- Using %zu bytes of GPU global memory.\n"
#define STRFMT_SETTINGS "- Running %u islands with %u individuals each.\n"
#define ERR_NOMEM "Couldn't allocate the required GPU memory.\n"
#define ERR_COSTKERN "Error in cost matrix kernel.\n"
#define ERR_INITPOP "Error while initializing the starting population.\n"
#define OK_CUTOFF "Stopped because we reached the cutoff point :)\n"
#define OK_EFFORT "Stopped because we reached the maximum effort\n"

static_assert(sizeof(v2) == sizeof(float2), "The size of float2 /= size of v2");

#define MUTATION_CHANCE 0.2f

__device__ inline r32
d_Distance(float2 U, float2 V)
{
    r32 Dx = U.x - V.x;
    r32 Dy = U.y - V.y;
    return __fsqrt_rn(fmaf(Dx, Dx, Dy * Dy));
}

__global__ void
FillCostMatrix(i32 N, float2 *Coords, r32 *CostMtx, i32 TileSize)
{
    /*
     * NOTE(srp): This could be optimized so that tiles are longer and each
     *            thread loads more than one coordinate into it, maximizing
     *            shared memory occupancy. At the moment it's overkill.
     */
    extern __shared__ float2 shm_Tile[];
    i32 T = threadIdx.x;
    i32 B = blockIdx.x;
    i32 I = T + blockDim.x*B;
    float2 RowCoord = Coords[I];

    i32 TileCount = (N + TileSize - 1)/TileSize;
    for (i32 TileIdx = 0; TileIdx < TileCount; TileIdx++)
    {
        i32 StartCol = TileSize*TileIdx;
        if (StartCol+T < N)
        {
            shm_Tile[T] = Coords[StartCol+T];
        }
        else
        {
            shm_Tile[T] = {};
        }
        __syncthreads();

        for (i32 J = 0; J < TileSize; J++)
        {
            i32 Col = StartCol+J;
            if (Col >= N) break;

            float2 ColCoord = shm_Tile[J];
            CostMtx[Col + N*I] = d_Distance(RowCoord, ColCoord);
        }
        __syncthreads();
    }
}

__global__ void
InitPopulation(i32 N, r32 *Costs, u32 *Islands, u32 *Elites, r32 *Fitness,
        r32 *EliteFitness, curandState *Rand, u64 Seed)
{
    // NOTE(srp): We need at least 6*blockDim.x bytes of shared memory per block
    extern __shared__ r32 Shm[];
    i32 Population = blockDim.x;
    r32 *shm_Fitness = Shm;
    u16 *shm_Indices = (u16*)(shm_Fitness + Population);

    i32 T = threadIdx.x;
    i32 B = blockIdx.x;
    i32 I = T + blockDim.x * B;
    u32 *Entry = Islands + I*N;

    curand_init(Seed + I, I, 0, &Rand[I]);
    curandState Rng = Rand[I];

    for (i32 K = 0; K < N; K++)
    {
        Entry[K] = K;
    }

    for (i32 Right = N-1; Right > 0; Right--)
    {
        i32 Left = curand(&Rng) % (Right + 1);
        i32 Tmp = Entry[Left];
        Entry[Left] = Entry[Right];
        Entry[Right] = Tmp;
    }

    r32 Evaluation = 0.0f;
    for (i32 K = 0; K < N; K++)
    {
        Evaluation += Costs[Entry[K] + N*Entry[(K+1)%N]];
    }
    shm_Fitness[T] = Evaluation;
    shm_Indices[T] = T;
    Fitness[I] = Evaluation;
    __syncthreads();

    for (i32 Stride = Population >> 1; Stride > 0; Stride >>= 1)
    {
        if (T < Stride)
        {
            if (shm_Fitness[T] > shm_Fitness[T+Stride])
            {
                shm_Fitness[T] = shm_Fitness[T+Stride];
                shm_Indices[T] = shm_Indices[T+Stride];
            }
        }
        __syncthreads();
    }

    __syncthreads();
    u32 *Best = Islands + N*(shm_Indices[0] + B*Population);
    for (i32 J = T; J < N; J += Population)
    {
        Elites[J + B*N] = Best[J];
    }

    Rand[I] = Rng;

    if (T == shm_Indices[0])
    {
        EliteFitness[B] = Evaluation;
    }
}

__global__ void
NextGeneration(i32 N, r32 *Costs, u32 *Islands, u32 *NewIslands, u32 *Elites,
        r32 *Fitness, r32 *EliteFitness, curandState *Rand, u64 *Bitsets, u32 BitsetGroups)
{
    // NOTE(srp): We need at least 6*blockDim.x bytes of shared memory per block
    extern __shared__ r32 Shm[];
    i32 Population = blockDim.x;
    r32 *shm_Fitness = Shm;
    u16 *shm_Indices = (u16*)(shm_Fitness + Population);

    i32 T = threadIdx.x;
    i32 B = blockIdx.x;
    i32 I = T + Population * B;

    u32 *OldEntry = Islands + I*N;
    u32 *NextEntry = NewIslands + I*N;

    curandState Rng = Rand[I];
    shm_Fitness[T] = Fitness[I];
    shm_Indices[T] = T;

    __syncthreads();

    i32 IslandCount = gridDim.x;
    i32 Candidate = curand(&Rng) % (IslandCount + Population);
    r32 CandidateFit = 0.0f;
    if (Candidate < Population)
    {
        CandidateFit = shm_Fitness[Candidate];
    }
    else
    {
        CandidateFit = EliteFitness[Candidate-Population];
    }

    for (i16 J = 0; J < 4; J++)
    {
        i32 NewCandidate = curand(&Rng) % (IslandCount + Population);
        if (NewCandidate < Population && shm_Fitness[NewCandidate] < shm_Fitness[Candidate])
        {
            Candidate = NewCandidate;
            CandidateFit = shm_Fitness[Candidate];
        }
        else
        {
            r32 NewFit = EliteFitness[NewCandidate-Population];
            if (NewFit < CandidateFit)
            {
                Candidate = NewCandidate;
                CandidateFit = NewFit;
            }
        }
    }

    u32 *Parent1 = Candidate < Population
        ? Islands + (Candidate + Population*B)*N
        : Elites + (Candidate-Population)*N;
    u32 *Parent2 = OldEntry;

    u32 P1 = curand(&Rng) % N; 
    u32 P2 = curand(&Rng) % N; 
    u32 CrossStart = min(P1, P2);
    u32 CrossEnd = max(P1, P2);

    u64 *MyBitset = Bitsets + I*BitsetGroups;
    for (i32 J = 0; J < BitsetGroups; J++)
    {
         MyBitset[J] = 0;
    }

    for (u32 J = CrossStart; J < CrossEnd; J++)
    {
         u32 K = Parent1[J];
         MyBitset[K/64] |= 1ull << (K%64);
    }

    u32 CrossIdx = 0;
    for (i32 J = 0; J < CrossEnd; J++)
    {
        u32 K = Parent2[J];
        if ((MyBitset[(K/64)] >> (K%64)) & 1)
        {
            continue;
        }
        NextEntry[CrossIdx++] = K;
    }
    for (i32 J = CrossStart; J < CrossEnd; J++)
    {
        NextEntry[CrossIdx++] = Parent1[J];
    }
    for (i32 J = CrossEnd; J < N; J++)
    {
        u32 K = Parent2[J];
        if ((MyBitset[(K/64)] >> (K%64)) & 1)
        {
            continue;
        }
        NextEntry[CrossIdx++] = K;
    }

    if (curand_uniform(&Rng) <= MUTATION_CHANCE)
    {
        P1 = curand(&Rng) % N; 
        P2 = curand(&Rng) % N; 
        CrossStart = min(P1, P2);
        CrossEnd = max(P1, P2);

        for (; CrossStart < CrossEnd; CrossStart++, CrossEnd--)
        {
            r32 Tmp = NextEntry[CrossStart];
            NextEntry[CrossStart] = NextEntry[CrossEnd];
            NextEntry[CrossEnd] = Tmp;
        }
    }

    r32 Evaluation = 0.0f;
    for (i32 K = 0; K < N; K++)
    {
        Evaluation += Costs[NextEntry[K] + N*NextEntry[(K+1)%N]];
    }

    if (Evaluation < shm_Fitness[T])
    {
        for (u32 K = 0; K < N; K++)
        {
            OldEntry[K] = NextEntry[K];
        }
        shm_Fitness[T] = Evaluation;
        Fitness[I] = Evaluation;
    }

    __syncthreads();

    for (i32 Stride = Population >> 1; Stride > 0; Stride >>= 1)
    {
        if (T < Stride)
        {
            if (shm_Fitness[T] > shm_Fitness[T+Stride])
            {
                shm_Fitness[T] = shm_Fitness[T+Stride];
                shm_Indices[T] = shm_Indices[T+Stride];
            }
        }
        __syncthreads();
    }

    __syncthreads();
    if (shm_Fitness[0] < EliteFitness[B])
    {
        u32 *Best = Islands + N*(shm_Indices[0] + B*Population);
        for (i32 J = T; J < N; J += Population)
        {
            Elites[J + B*N] = Best[J];
        }

        if (T == shm_Indices[0])
        {
            EliteFitness[B] = Evaluation;
        }
    }

    Rand[I] = Rng;
}

internal inline b32
Main(tsp_instance *__restrict__ Tsp,
        i32 *__restrict__ out_Permutation,
        u64 *__restrict__ Iterations,
        cudaDeviceProp *__restrict__ DevProp,
        r32 Cutoff)
{
    const auto SmCount = DevProp->multiProcessorCount;
    const auto MaxSmThreads = DevProp->maxThreadsPerMultiProcessor;
    const auto MaxBlockThreads = DevProp->maxThreadsPerBlock;
    const auto ShmPerBlock = DevProp->sharedMemPerBlock;

    // The strategy is to fill up the entire GPU's compute availability exactly
    // a fraction of N times each kernel call to reduce kernel launch
    // synchronization overhead.
    const u32 Population = (max(
            min(MaxSmThreads / SmCount, MaxBlockThreads),
            4*DevProp->warpSize)/DevProp->warpSize)*DevProp->warpSize;
    const u32 Islands = SmCount * ((Tsp->N+Population-1)/Population);
    char SettingsStr[sizeof(STRFMT_SETTINGS) + 21] = {};
    size_t SettingsStrLen = sprintf(SettingsStr, STRFMT_SETTINGS, Islands, Population);
    write(1, SettingsStr, SettingsStrLen);

    void *d_Arena;
    r32 *d_CostMtx, *d_Fitness, *d_EliteFitness;
    float2 *d_Coords;
    u32 *d_IslandsA, *d_IslandsB, *d_EliteIsland;
    curandState *d_Rand;
    u64 *d_Bitsets;

    const size_t d_CostMtxSz = sizeof(r32) * Tsp->N * Tsp->N;
    const size_t d_FitnessSz = sizeof(r32) * Islands * Population;
    const size_t d_EliteFitnessSz = sizeof(r32) * Islands;
    const size_t d_CoordsSz = sizeof(float2) * Tsp->N;
    const size_t d_IslandsASz = sizeof(u32) * Islands * Population * Tsp->N;
    const size_t d_IslandsBSz = sizeof(u32) * Islands * Population * Tsp->N;
    const size_t d_EliteIslandSz = sizeof(u32) * Islands * Tsp->N;
    const size_t d_RandSz = sizeof(curandState) * Islands * Population;
    const size_t d_BitsetsSz = sizeof(u64) * Islands * Population * ((Tsp->N+63)/64);

    const size_t ArenaSize = d_CostMtxSz + d_FitnessSz + d_EliteFitnessSz
        + d_CoordsSz + d_IslandsASz + d_IslandsBSz + d_EliteIslandSz
        + d_RandSz + d_BitsetsSz;
    char MemUseStr[sizeof(STRFMT_GLOBUSE) + 21] = {};
    size_t MemUseStrLen = sprintf(MemUseStr, STRFMT_GLOBUSE, ArenaSize);
    write(1, MemUseStr, MemUseStrLen);

    if (ArenaSize > DevProp->totalGlobalMem)
    {
        write(1, ERR_NOMEM, sizeof(ERR_NOMEM));
        return 0;
    }

    cudaError_t Err = cudaMalloc(&d_Arena, ArenaSize);
    if (Err != cudaSuccess)
    {
        write(1, ERR_NOMEM, sizeof(ERR_NOMEM));
        return 0;
    }

    d_Bitsets = (u64*)d_Arena;
    d_Rand = (curandState*)(((u8*)d_Bitsets)+d_BitsetsSz);
    d_Coords = (float2*)(((u8*)d_Rand)+d_RandSz);
    d_CostMtx = (r32*)(((u8*)d_Coords)+d_CoordsSz);
    d_IslandsA = (u32*)(((u8*)d_CostMtx)+d_CostMtxSz);
    d_IslandsB = (u32*)(((u8*)d_IslandsA)+d_IslandsASz);
    d_EliteIsland = (u32*)(((u8*)d_IslandsB)+d_IslandsBSz);
    d_Fitness = (r32*)(((u8*)d_EliteIsland)+d_EliteIslandSz);
    d_EliteFitness = (r32*)(((u8*)d_Fitness)+d_FitnessSz);

    Err = cudaMemcpy(d_Coords, Tsp->Coords, d_CoordsSz, cudaMemcpyHostToDevice);
    if (Err != cudaSuccess)
    {
        cudaFree(d_Arena);
        return 0;
    }

    i32 CostMtxMaxTileSize = ShmPerBlock/(sizeof(float2));
    i32 CostThreads = CostMtxMaxTileSize < MaxBlockThreads
        ? CostMtxMaxTileSize
        : MaxBlockThreads;
    i32 CostBlocks = (Tsp->N + CostThreads - 1)/CostThreads;
    i32 CostShm = sizeof(float2) * CostThreads;
    FillCostMatrix<<<CostBlocks,CostThreads,CostShm>>>(Tsp->N, d_Coords, d_CostMtx, CostThreads);

    cudaDeviceSynchronize();
    if ((Err = cudaGetLastError()) != cudaSuccess)
    {
        write(1, ERR_COSTKERN, sizeof(ERR_COSTKERN));
        cudaFree(d_Arena);
        return 0;
    }

    u32 ShmSize = 6*Population;
    InitPopulation<<<Islands,Population,ShmSize>>>(Tsp->N, d_CostMtx,
            d_IslandsA, d_EliteIsland, d_Fitness, d_EliteFitness, d_Rand,
            __rdtsc());

    cudaDeviceSynchronize();
    if ((Err = cudaGetLastError()) != cudaSuccess)
    {
        write(1, ERR_INITPOP, sizeof(ERR_INITPOP));
        cudaFree(d_Arena);
        return 0;
    }

    r32 h_EliteFitness[Islands];
    r32 BestFit = INFINITY;

    u64 I = 0;
    u64 Its = *Iterations;
    for (; I < Its && BestFit > Cutoff; I++)
    {
        NextGeneration<<<Islands, Population, ShmSize>>>(Tsp->N, d_CostMtx,
                d_IslandsA, d_IslandsB, d_EliteIsland, d_Fitness, d_EliteFitness,
                d_Rand, d_Bitsets, (Tsp->N+63)/64);
        cudaMemcpy(h_EliteFitness, d_EliteFitness, sizeof(h_EliteFitness), cudaMemcpyDeviceToHost);
        i32 BestJ = -1;
        for (u32 J = 0; J < Islands; J++)
        {
            if (h_EliteFitness[J] < BestFit)
            {
                BestJ = J;
                BestFit = h_EliteFitness[J];
            }
        }
        if (BestJ >= 0)
        {
            cudaMemcpy(out_Permutation, d_EliteIsland + BestJ*Tsp->N,
                    sizeof(u32)*Tsp->N, cudaMemcpyDeviceToHost);
        }
    }

    if (BestFit <= Cutoff)
    {
        write(1, OK_CUTOFF, sizeof(OK_CUTOFF));
    }
    else
    {
        write(1, OK_EFFORT, sizeof(OK_EFFORT));
    }

    *Iterations = I;

    cudaFree(d_Arena);
    return 1;
}

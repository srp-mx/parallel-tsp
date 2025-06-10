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

/**
 * Device function that calculates the euclidean distance quickly.
 *
 * @param U The first coordinate.
 * @param V The second coordinate.
 *
 * @return The distance from U to V in the L2 norm.
 */
__device__ inline r32
d_Distance(float2 U, float2 V)
{
    r32 Dx = U.x - V.x;
    r32 Dy = U.y - V.y;
    return __fsqrt_rn(fmaf(Dx, Dx, Dy * Dy));
}

/**
 * CUDA kernel that fills in the cost matrix.
 *
 * @param N The number of coordinates.
 * @param Coords Array of 2D coordinates.
 * @param CostMtx Space of r32's to fill in with the matrix's entries.
 * @param TileSize Byte size of the shared memory tiles.
 */
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

/**
 * CUDA kernel that initializes the starting population.
 *
 * @param N The number of cities.
 * @param Costs The NxN cost matrix in row-major order.
 * @param Islands Pointer to the population individuals.
 * @param Elites Pointer to the elite population individuals.
 * @param Fitness The fitness of each individual.
 * @param EliteFitness The fitness of each elite individual.
 * @param Rand Array of random generators.
 * @param Seed Random generation base seed.
 */
__global__ void
InitPopulation(i32 N,
        r32 *Costs,
        u32 *Islands,
        u32 *Elites,
        r32 *Fitness,
        r32 *EliteFitness,
        curandState *Rand,
        u64 Seed)
{
    // NOTE(srp): We need at least 6*blockDim.x bytes of shared memory per block
    extern __shared__ r32 Shm[];
    i32 Population = blockDim.x;
    r32 *shm_Fitness = Shm;
    u16 *shm_Indices = (u16*)(shm_Fitness + Population);

    i32 T = threadIdx.x;
    i32 B = blockIdx.x;
    i32 I = T + blockDim.x * B;

    // Points to this thread's individual
    u32 *Entry = Islands + I*N;

    // Initializes the random state of each generator (one per individual)
    curand_init(Seed + I, I, 0, &Rand[I]);
    curandState Rng = Rand[I];

    // Start the permutation by filling the entries from 0 to N-1
    for (i32 K = 0; K < N; K++)
    {
        Entry[K] = K;
    }

    // Shuffle the entries to get a valid permutation
    for (i32 Right = N-1; Right > 0; Right--)
    {
        i32 Left = curand(&Rng) % (Right + 1);
        i32 Tmp = Entry[Left];
        Entry[Left] = Entry[Right];
        Entry[Right] = Tmp;
    }

    // Evaluate and share with the thread block
    r32 Evaluation = 0.0f;
    for (i32 K = 0; K < N; K++)
    {
        Evaluation += Costs[Entry[K] + N*Entry[(K+1)%N]];
    }
    shm_Fitness[T] = Evaluation;
    shm_Indices[T] = T;
    Fitness[I] = Evaluation;
    __syncthreads();

    // Find optimal individual of the thread block through parallel reduction
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

    // Store the best individual of the thread block on the elite island
    __syncthreads();
    u32 *Best = Islands + N*(shm_Indices[0] + B*Population);
    for (i32 J = T; J < N; J += Population)
    {
        Elites[J + B*N] = Best[J];
    }

    // Write back the RNG state to avoid repetition
    Rand[I] = Rng;

    // Let the winning thread write back its fitness from its registers
    if (T == shm_Indices[0])
    {
        EliteFitness[B] = Evaluation;
    }
}

/**
 * CUDA kernel that builds the next generation from the current one.
 *
 * @param N The number of cities.
 * @param Costs The cost matrix in row-major order.
 * @param Islands Pointer to the population individuals, to be replaced.
 * @param NewIslands Pointer to a buffer with all new candidate individuals.
 * @param Elites Pointer to the elite population individuals.
 * @param Fitness The fitness of each individual.
 * @param EliteFitness The fitness of each elite individual.
 * @param Rand Array of random generators.
 * @param Bitsets Buffer per individual with some u64 to use as bitsets when
 *                doing the crossover operation.
 * @param BitsetGroups How many u64 are given to each thread.
 */
__global__ void
NextGeneration(i32 N,
        r32 *Costs,
        u32 *Islands,
        u32 *NewIslands,
        u32 *Elites,
        r32 *Fitness,
        r32 *EliteFitness,
        curandState *Rand,
        u64 *Bitsets,
        u32 BitsetGroups)
{
    // NOTE(srp): We need at least 6*blockDim.x bytes of shared memory per block
    extern __shared__ r32 Shm[];
    i32 Population = blockDim.x;
    r32 *shm_Fitness = Shm;
    u16 *shm_Indices = (u16*)(shm_Fitness + Population);

    i32 T = threadIdx.x;
    i32 B = blockIdx.x;
    i32 I = T + Population * B;

    // Points to the thread's individual and buffer for the new one
    u32 *OldEntry = Islands + I*N;
    u32 *NextEntry = NewIslands + I*N;

    // Fetches the cuRand state from global memory
    curandState Rng = Rand[I];
    shm_Fitness[T] = Fitness[I];
    shm_Indices[T] = T;

    __syncthreads();

    // Gets the first candidate for the tournament selection operation
    i32 IslandCount = gridDim.x;
    i32 Candidate = curand(&Rng) % (IslandCount + Population);
    r32 CandidateFit = 0.0f;
    if (Candidate < Population)
    {
        // We take from our population
        CandidateFit = shm_Fitness[Candidate];
    }
    else
    {
        // We take from the elite island
        CandidateFit = EliteFitness[Candidate-Population];
    }

    // Gets the remaining 4 candidates for tournament selection, and keeps only
    // the best candidate seen so far.
    for (i16 J = 0; J < 4; J++)
    {
        i32 NewCandidate = curand(&Rng) % (IslandCount + Population);
        if (NewCandidate < Population && shm_Fitness[NewCandidate] < shm_Fitness[Candidate])
        {
            // We take from our population
            Candidate = NewCandidate;
            CandidateFit = shm_Fitness[Candidate];
        }
        else
        {
            r32 NewFit = EliteFitness[NewCandidate-Population];
            if (NewFit < CandidateFit)
            {
                // We take from the elite island
                Candidate = NewCandidate;
                CandidateFit = NewFit;
            }
        }
    }

    // We get the pointer for the parent chosen, whether from this island or
    // the elite island.
    u32 *Parent1 = Candidate < Population
        ? Islands + (Candidate + Population*B)*N
        : Elites + (Candidate-Population)*N;
    // The individual who lives in this thread will be the other parent.
    u32 *Parent2 = OldEntry;

    // We get the two-point crossover indices
    u32 P1 = curand(&Rng) % N; 
    u32 P2 = curand(&Rng) % N; 
    u32 CrossStart = min(P1, P2);
    u32 CrossEnd = max(P1, P2);

    // We get a pointer to this thread's bitset and clear it out
    u64 *MyBitset = Bitsets + I*BitsetGroups;
    for (i32 J = 0; J < BitsetGroups; J++)
    {
         MyBitset[J] = 0;
    }

    // Load the crossover elements from the middle parent into the bitset
    for (u32 J = CrossStart; J < CrossEnd; J++)
    {
         u32 K = Parent1[J];
         MyBitset[K/64] |= 1ull << (K%64);
    }

    // Compute the crossover operation
    u32 CrossIdx = 0;
    for (i32 J = 0; J < CrossEnd; J++)
    {
        u32 K = Parent2[J];
        // We add from the second (outside) parent only if the current element
        // is not in the bitset. We add from 0 up to the crossover end index
        // since we must add every city.
        if ((MyBitset[(K/64)] >> (K%64)) & 1)
        {
            continue;
        }
        NextEntry[CrossIdx++] = K;
    }
    for (i32 J = CrossStart; J < CrossEnd; J++)
    {
        // We add all elements from the first (inside) parent.
        NextEntry[CrossIdx++] = Parent1[J];
    }
    for (i32 J = CrossEnd; J < N; J++)
    {
        u32 K = Parent2[J];
        // We add from the second (outside) parent only if the current element
        // is not in the bitset. We add from the crossover end index to the end
        // since we must add every city.
        if ((MyBitset[(K/64)] >> (K%64)) & 1)
        {
            continue;
        }
        NextEntry[CrossIdx++] = K;
    }

    // We get a chance to mutate randomly
    if (curand_uniform(&Rng) <= MUTATION_CHANCE)
    {
        // Chose two indices to mutate
        P1 = curand(&Rng) % N; 
        P2 = curand(&Rng) % N; 
        CrossStart = min(P1, P2);
        CrossEnd = max(P1, P2);

        // Invert the range defined by those two indices, that makes the
        // mutation only change 2 adjacencies.
        for (; CrossStart < CrossEnd; CrossStart++, CrossEnd--)
        {
            r32 Tmp = NextEntry[CrossStart];
            NextEntry[CrossStart] = NextEntry[CrossEnd];
            NextEntry[CrossEnd] = Tmp;
        }
    }

    // Make each thread evaluate its individual
    r32 Evaluation = 0.0f;
    for (i32 K = 0; K < N; K++)
    {
        Evaluation += Costs[NextEntry[K] + N*NextEntry[(K+1)%N]];
    }

    // If the new individual is better, or it passes a random check, it
    // replaces the old individual in the population.
    if (Evaluation < shm_Fitness[T] || curand_uniform(&Rng) <= ACCEPT_WORSE_CHANCE)
    {
        for (u32 K = 0; K < N; K++)
        {
            OldEntry[K] = NextEntry[K];
        }
        shm_Fitness[T] = Evaluation;
        Fitness[I] = Evaluation;
    }

    __syncthreads();

    // We get the best individual from this island by way of parallel reduction
    for (i32 Stride = Population >> 1; Stride > 0; Stride >>= 1)
    {
        if (T < Stride && shm_Fitness[T] > shm_Fitness[T+Stride])
        {
            shm_Fitness[T] = shm_Fitness[T+Stride];
            shm_Indices[T] = shm_Indices[T+Stride];
        }
        __syncthreads();
    }

    // Have thread 0 from this block load this island's old elite's fitness
    if (T == 0)
    {
        shm_Fitness[1] = EliteFitness[B];
    }
    __syncthreads();
    // If this elite's fitness is better than the last, we replace it in
    // parallel.
    if (shm_Fitness[0] < shm_Fitness[1])
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

    // Write back the random state to avoid repetition
    Rand[I] = Rng;
}

/**
 * Utility structure to pass as a unary function for the CUDA API.
 */
struct ConvertBlockSizeToShmSize
{
    /**
     * Returns the amount of shared memory needed given a particular thread
     * block size. This function is intended for the NextGeneration and perhaps
     * InitPopulation kernels, and it aids CUDA estimate ideal block and grid
     * size to maximize occupancy.
     *
     * @param BlockSize The number of threads in a block.
     *
     * @return The number of bytes of shared memmory required given the block
     *         size.
     */
    __host__ __device__ size_t
    operator()(i32 BlockSize)
    {
        return 6 * BlockSize;
    }
};

/**
 * Runs the solver.
 *
 * @param Tsp Copy of the problem instance.
 * @param out_Permutation Pointer to the array which will hold the result.
 * @param Iterations Pointer to the number of iterations, which must be written
 *                   back to if the actual iterations were less.
 * @param DevProp Device properties object.
 * @param Cutoff If the costs reach this or lower, we may end early.
 *
 * @return 1 if everything went ok, 0 otherwise.
 */
internal inline b32
Main(tsp_instance *__restrict__ Tsp,
        i32 *__restrict__ out_Permutation,
        u64 *__restrict__ Iterations,
        cudaDeviceProp *__restrict__ DevProp,
        r32 Cutoff,
        i32 Parallelism)
{
    const auto MaxBlockThreads = DevProp->maxThreadsPerBlock;
    const auto ShmPerBlock = DevProp->sharedMemPerBlock;
    const auto WarpSize = DevProp->warpSize;

    // Maximize occupancy for a net population size of at least 2N
    u32 Population = 0;
    u32 Islands = 0;
    cudaOccupancyMaxPotentialBlockSizeVariableSMem(
            (i32*)&Islands, (i32*)&Population, NextGeneration,
            ConvertBlockSizeToShmSize());

    if (Parallelism == 0)
    {
        Population = max(4*WarpSize, Population);
        Islands = max(Islands, (2*Tsp->N + Population - 1)/Population);
    }
    else
    {
        Population = max(4*WarpSize, Population);
        Islands = Parallelism;
    }

    // Print execution settings
    char SettingsStr[sizeof(STRFMT_SETTINGS) + 21] = {};
    size_t SettingsStrLen = sprintf(SettingsStr, STRFMT_SETTINGS, Islands, Population);
    IGNORE_RESULT(write(1, SettingsStr, SettingsStrLen));

    // Pointers to the different buffers in global device memory
    void *d_Arena;
    r32 *d_CostMtx, *d_Fitness, *d_EliteFitness;
    float2 *d_Coords;
    u32 *d_IslandsA, *d_IslandsB, *d_EliteIsland;
    curandState *d_Rand;
    u64 *d_Bitsets;

    // The sizes of each of these buffers
    const size_t d_CostMtxSz = sizeof(r32) * Tsp->N * Tsp->N;
    const size_t d_FitnessSz = sizeof(r32) * Islands * Population;
    const size_t d_EliteFitnessSz = sizeof(r32) * Islands;
    const size_t d_CoordsSz = sizeof(float2) * Tsp->N;
    const size_t d_IslandsASz = sizeof(u32) * Islands * Population * Tsp->N;
    const size_t d_IslandsBSz = sizeof(u32) * Islands * Population * Tsp->N;
    const size_t d_EliteIslandSz = sizeof(u32) * Islands * Tsp->N;
    const size_t d_RandSz = sizeof(curandState) * Islands * Population;
    const size_t d_BitsetsSz = sizeof(u64) * Islands * Population * ((Tsp->N+63)/64);

    // Size of the net allocation, to be done at once on a single arena.
    // Note that the drawbacks of this is that if memory consumption is high,
    // we may fail even if there is enough available memory on grounds of
    // fragmentation. However, this makes resource handling a little less
    // cumbersome as well as possibly more efficient due to the contiguity.
    const size_t ArenaSize = d_CostMtxSz + d_FitnessSz + d_EliteFitnessSz
        + d_CoordsSz + d_IslandsASz + d_IslandsBSz + d_EliteIslandSz
        + d_RandSz + d_BitsetsSz;
    char MemUseStr[sizeof(STRFMT_GLOBUSE) + 21] = {};
    size_t MemUseStrLen = sprintf(MemUseStr, STRFMT_GLOBUSE, ArenaSize);
    IGNORE_RESULT(write(1, MemUseStr, MemUseStrLen));

    // If we need more memory than available, fail and tell the user why
    if (ArenaSize > DevProp->totalGlobalMem)
    {
        IGNORE_RESULT(write(1, ERR_NOMEM, sizeof(ERR_NOMEM)));
        return 0;
    }

    // If we couldn't allocate the arena, fail and tell the user why
    cudaError_t Err = cudaMalloc(&d_Arena, ArenaSize);
    if (Err != cudaSuccess)
    {
        IGNORE_RESULT(write(1, ERR_NOMEM, sizeof(ERR_NOMEM)));
        return 0;
    }

    // Split the arena up
    d_Bitsets = (u64*)d_Arena;
    d_Rand = (curandState*)(((u8*)d_Bitsets)+d_BitsetsSz);
    d_Coords = (float2*)(((u8*)d_Rand)+d_RandSz);
    d_CostMtx = (r32*)(((u8*)d_Coords)+d_CoordsSz);
    d_IslandsA = (u32*)(((u8*)d_CostMtx)+d_CostMtxSz);
    d_IslandsB = (u32*)(((u8*)d_IslandsA)+d_IslandsASz);
    d_EliteIsland = (u32*)(((u8*)d_IslandsB)+d_IslandsBSz);
    d_Fitness = (r32*)(((u8*)d_EliteIsland)+d_EliteIslandSz);
    d_EliteFitness = (r32*)(((u8*)d_Fitness)+d_FitnessSz);

    // If memcpy somehow fails, just fail since that would be odd
    Err = cudaMemcpy(d_Coords, Tsp->Coords, d_CoordsSz, cudaMemcpyHostToDevice);
    if (Err != cudaSuccess)
    {
        cudaFree(d_Arena);
        return 0;
    }

    // Calculate parameters for the cost matrix construction kernel
    i32 CostMtxMaxTileSize = ShmPerBlock/(sizeof(float2));
    i32 CostThreads = CostMtxMaxTileSize < MaxBlockThreads
        ? CostMtxMaxTileSize
        : MaxBlockThreads;
    i32 CostBlocks = (Tsp->N + CostThreads - 1)/CostThreads;
    i32 CostShm = sizeof(float2) * CostThreads;
    // Fill the cost matrix
    FillCostMatrix<<<CostBlocks,CostThreads,CostShm>>>(
            Tsp->N,
            d_Coords,
            d_CostMtx,
            CostThreads);

    // Check if the computation failed, in which case fail and tell the user
    cudaDeviceSynchronize();
    if ((Err = cudaGetLastError()) != cudaSuccess)
    {
        IGNORE_RESULT(write(1, ERR_COSTKERN, sizeof(ERR_COSTKERN)));
        cudaFree(d_Arena);
        return 0;
    }

    // Build the initial population
    u32 ShmSize = ConvertBlockSizeToShmSize()(Population);
    u64 RandSeed = __rdtsc();
    InitPopulation<<<Islands,Population,ShmSize>>>(
            Tsp->N,
            d_CostMtx,
            d_IslandsA,
            d_EliteIsland,
            d_Fitness,
            d_EliteFitness,
            d_Rand,
            RandSeed);

    // Check if the computation failed, in which case fail and tell the user
    cudaDeviceSynchronize();
    if ((Err = cudaGetLastError()) != cudaSuccess)
    {
        IGNORE_RESULT(write(1, ERR_INITPOP, sizeof(ERR_INITPOP)));
        cudaFree(d_Arena);
        return 0;
    }

    // Set up a buffer to copy elite fitnesses on each iteration in order to
    // compare them and find the best fitness and individual.
    r32 h_EliteFitness[Islands];
    r32 BestFit = INFINITY;

    u64 I = 0;
    u64 Its = *Iterations;
    for (; I < Its && BestFit > Cutoff; I++)
    {
        // Build the next generation
        NextGeneration<<<Islands, Population, ShmSize>>>(
                Tsp->N,
                d_CostMtx,
                d_IslandsA,
                d_IslandsB,
                d_EliteIsland,
                d_Fitness,
                d_EliteFitness,
                d_Rand,
                d_Bitsets,
                (Tsp->N+63)/64);
        // Copy the elite island into host memory
        cudaMemcpy(h_EliteFitness, d_EliteFitness, sizeof(h_EliteFitness),
                cudaMemcpyDeviceToHost);
        // Check if the best solution found is here
        i32 BestJ = -1;
        for (u32 J = 0; J < Islands; J++)
        {
            if (h_EliteFitness[J] < BestFit)
            {
                BestJ = J;
                BestFit = h_EliteFitness[J];
            }
        }

        // If we found our new best solution, copy the route into the result
        // buffer (out_Permutation).
        if (BestJ >= 0)
        {
            cudaMemcpy(out_Permutation, d_EliteIsland + BestJ*Tsp->N,
                    sizeof(u32)*Tsp->N, cudaMemcpyDeviceToHost);
        }
    }

    if (BestFit <= Cutoff)
    {
        // If the best solution was under or at the cutoff, we tell the user
        IGNORE_RESULT(write(1, OK_CUTOFF, sizeof(OK_CUTOFF)));
    }
    else
    {
        // Otherwise, we tell the user it didn't
        IGNORE_RESULT(write(1, OK_EFFORT, sizeof(OK_EFFORT)));
    }

    // Write back the number of iterations executed
    *Iterations = I;

    // Free device memory
    cudaFree(d_Arena);

    // Exit successfully
    return 1;
}

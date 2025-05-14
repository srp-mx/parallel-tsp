#pragma once

#include "defines.h"
#include "tsp.h"

#include <omp.h>

#include <unistd.h>
#include <stdio.h>
#include <sys/mman.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <random>
#include <immintrin.h>
#include <algorithm>

#define ERR_MALLOC "No se pudo alojar suficiente memoria.\n"
#define NUM_ELITE 5
#define MUTACION 0.1f
#define PROBA_ELITISMO 0.05f
#define PROBA_HIJO_GRATIS 0.1f

struct rng
{
    std::random_device RD_;
    std::mt19937 Generator_;

    rng()
      : Generator_(RD_())
    { }

    i32
    Rng_i32(i32 N)
    {
        std::uniform_int_distribution<i32> Distribution(0, N-1);
        return Distribution(Generator_);
    }

    r32
    Rng_r32()
    {
        std::uniform_real_distribution<r32> Distribution(0.0f, 1.0f);
        return Distribution(Generator_);
    }
};

/**
 * Obtiene el apuntador a un individuo de la población.
 *
 * @param N El número de ciudades.
 * @param I El índice del individuo.
 *
 * @return El apuntador a la primera ciudad del individuo.
 */
internal inline i32*
GetIndividuo(i32 *Poblacion, i32 N, i32 I)
{
    return Poblacion + I*N;
}

/**
 * Obtiene el apuntador a una población de una lista de islas.
 *
 * @param P El número de individuos en cada población/isla.
 * @param N El número de ciudades.
 * @param I El índice del individuo.
 *
 * @return El apuntador a la primera ciudad del individuo.
 */
internal inline i32*
GetPoblacion(i32 *Islas, i32 P, i32 N, i32 I)
{
    return Islas + P*N*I;
}

template <typename T>
internal inline T*
GetDatoPoblacion(T *Datos, i32 P, i32 I)
{
    return Datos + P*I;
}

/**
 * Crea una población aleatoria de P soluciones de tamaño N.
 *
 * @param Poblacion El buffer con los N*P i32 a llenar.
 * @param N Tamaño de las soluciones o individuos.
 * @param P Número de individuos de la población.
 */
internal void
IniciaPoblacion(i32 *Poblacion, i32 N, i32 P, rng *Rng)
{
    for (i32 I = 0; I < P; I++)
    {
        i32 *Individuo = GetIndividuo(Poblacion, N, I);
        // Llenamos el arreglo con los números del 0 a N-1.
        for (i32 J = 0; J < N; J++)
        {
            Individuo[J] = J;
        }
        // Hacemos shuffle de los números.
        for (i32 Right = N-1; Right > 0; Right--)
        {
            i32 Left = Rng->Rng_i32(Right+1);
            i32 Tmp = Individuo[Left];
            Individuo[Left] = Individuo[Right];
            Individuo[Right] = Tmp;
        }
    }
}

/**
 * Calcula la distancia entre Ciudad1 y Ciudad2.
 *
 * @param Ciudad1 Primera ciudad.
 * @param Ciudad2 Segunda ciudad.
 * @param Coords Arreglo que contiene las coordenadas de todas las ciudades.
 * 
 * @return La distancia euclidiana entre las ciudades.
 */
internal r32
Distancia(i32 Ciudad1, i32 Ciudad2, v2 Coords[])
{
    r32 Dx = Coords[Ciudad1].X - Coords[Ciudad2].X;
    r32 Dy = Coords[Ciudad1].Y - Coords[Ciudad2].Y;
    return sqrtf(Dx*Dx + Dy*Dy);
}

/**
 * Calcula la aptitud de una permutación de ciudades.
 *
 * @param N Tamaño de las soluciones o individuos.
 * @param Individuo La solución de la que calculamos la aptitud.
 * @param Coords Arreglo que contiene las coordenadas de todas las ciudades.
 * 
 * @return Acc - La aptitud de una solución.
 */
internal r32
Aptitud(i32 N, i32 Individuo[], v2 Coords[])
{
    r32 Acc = 0.0f;
    for(i32 I = 0; I < N; I++)
    {
        Acc += Distancia(Individuo[I], Individuo[(I+1) % N], Coords);
    }
    return Acc;
}

/**
 * Consigue la posición en la que se encuentra la solución mejor evaluada.
 *
 * @param Puntuaciones Arreglo que contiene puntajes de las soluciones.
 * @param Len Número de elementos en el arreglo Puntuaciones.
 * 
 * @return La posición que contiene el menor valor en el arreglo Puntuaciones.
 */
internal i32
ConsigueMejor(r32 Puntuaciones[], i32 Len)
{
    r32 Menor = Puntuaciones[0];
    i32 Posicion = 0;
    for (i32 I = 1; I < Len; I++)
    {
        if (Puntuaciones[I] < Menor)
        {
            Menor = Puntuaciones[I];
            Posicion = I;
        }
    }
    return Posicion;
}

/**
 * Consigue el ganador de un torneo entre 5 soluciones aleatorias.
 * El ganador será el que tenga la calificación menor.
 *
 * @param Puntuaciones Arreglo que contiene puntajes de las soluciones.
 * @param Len Número de elementos en el arreglo Puntuaciones.
 * 
 * @return La posición que contiene el ganador del torneo.
 */
internal i32
Torneo(r32 Puntuaciones[], i32 Len, rng *Rng)
{
    i32 Random1 = Rng->Rng_i32(Len);
    i32 Random2 = Rng->Rng_i32(Len);
    i32 Random3 = Rng->Rng_i32(Len);
    i32 Random4 = Rng->Rng_i32(Len);
    i32 Random5 = Rng->Rng_i32(Len);
    i32 Pos[] = {Random1, Random2, Random3, Random4, Random5};
    r32 Arr[] = {
        Puntuaciones[Random1],
        Puntuaciones[Random2],
        Puntuaciones[Random3],
        Puntuaciones[Random4],
        Puntuaciones[Random5]
    };
    return Pos[ConsigueMejor(Arr, 5)];
}

/**
 * Toma dos soluciones padre y las combina para generar una nueva solución hija.
 *
 * @param NuevoIndividuo Donde se guardará el nuevo individuo que generamos.
 * @param Padre1 Solución del primer padre, la cual se recombinará.
 * @param Padre2 Solución del segundo padre, la cual se recombinará.
 * @param Mutacion Probabilidad de que la solución mute.
 * @param N Número de ciudades.
 */
internal void
Recombina(i32 *NuevoIndividuo, i32 *Padre1, i32 *Padre2, r32 Mutacion, i32 N, rng *Rng)
{
    //Combinamos a los padres.
    i32* Solucion = NuevoIndividuo;
    std::vector<bool> Bitmap(N, 0);
    i32 Pos1 = Rng->Rng_i32(N);
    i32 Pos2 = Rng->Rng_i32(N);
    i32 A = std::min(Pos1,Pos2);
    i32 B = std::max(Pos1,Pos2);
    for (i32 I = A; I < B; I++)
    {
        Bitmap[Padre1[I]] = 1;
    }
    i32 Idx = 0;
    for (i32 I = 0; I < B; I++)
    {
        if (Bitmap[Padre2[I]] == 0)
        {
            Solucion[Idx++] = Padre2[I];
        }
    }
    for (i32 I = A; I < B; I++)
    {
        Solucion[Idx++] = Padre1[I];
    }
    for (i32 I = B; I < N; I++)
    {
        if (Bitmap[Padre2[I]] == 0)
        {
            Solucion[Idx++] = Padre2[I];
        }
    }
    
    // Aplicamos mutación al invertir un intervalo.
    if (Rng->Rng_r32() <= Mutacion)
    {
        i32 Random1 = Rng->Rng_i32(N);
        i32 Random2 = Rng->Rng_i32(N);
        A = std::min(Random1, Random2);
        B = std::max(Random1, Random2);

        while (A < B)
        {
            i32 Temp = Solucion[A];
            Solucion[A] = Solucion[B];
            Solucion[B] = Temp;
            A++;
            B--;
        }
    }
}

/**
 * Encuentra los mejores NUM_ELITE valores entre la generación y las respuestas
 * previamente contenidas en la élite.
 *
 * @param Elite Un arreglo que contiene las NUM_ELITE mejores soluciones
 *              encontradas al momento.
 * @param EliteNueva Arreglo donde se guardará la nueva élite que generaremos.
 * @param Poblacion Arreglo con las soluciones de la generación actual.
 * @param Puntuaciones Puntos asignados a cada solución de la generación actual.
 */
internal void
ActualizaElite(i32 *Elite, i32 *EliteNueva, i32 *Poblacion, r32 *Puntuaciones,
        i32 P, i32 N, bool PrimeraVez)
{
    // Guardamos la puntuación y el índice de la solución en Poblacion.
    std::vector<std::pair<r32, i32>> Pares;
    for (i32 I = 0; I < P; ++I) {
        Pares.push_back({Puntuaciones[I], I});
    }

    // Ordenamos las puntuaciones en orden ascendente
    //(dando prioridad a las menores puntuaciones).
    std::sort(Pares.begin(), Pares.end());

    // Guardamos todos los individuos, tanto de la élite como de la población.
    std::vector<i32*> CandidatosAll;
    std::vector<r32> PuntuacionesAll;
    
    // Agregamos a la élite
    for (i32 I = 0; I < NUM_ELITE; ++I) {
        if (!PrimeraVez) {
            CandidatosAll.push_back(GetIndividuo(Elite, N, I));
            PuntuacionesAll.push_back(std::numeric_limits<r32>::max());
        }
    }

    // Agregamos la población normal
    i32 Candidatos = std::min(P, NUM_ELITE);
    for (i32 I = 0; I < Candidatos; ++I) {
        i32 Idx = Pares[I].second;
        CandidatosAll.push_back(GetIndividuo(Poblacion, N, Idx));
        PuntuacionesAll.push_back(Pares[I].first);
    }

    // Creamos un vector de pares para ordenar todos los candidatos.
    std::vector<std::pair<r32, i32*>> CandidatosOrdenados;
    for (size_t I = 0; I < CandidatosAll.size(); ++I) {
        CandidatosOrdenados.push_back({PuntuacionesAll[I], CandidatosAll[I]});
    }

    // Ordenamos a los candidatos.
    std::sort(CandidatosOrdenados.begin(), CandidatosOrdenados.end());

    // Copiamos las 5 mejores soluciones a la nueva élite.
    for (i32 I = 0; I < NUM_ELITE; ++I) {
        i32 *IndividuoViejo = CandidatosOrdenados[I].second;
        i32 *IndividuoNuevo = GetIndividuo(EliteNueva, N, I);
        memcpy(IndividuoNuevo, IndividuoViejo, sizeof(i32)*N);
    }
}

internal void
RemezclaElite(i32 N, i32 *Elite, r32 *PuntuacionesElite)
{
    i32 Buf[N];
    i32 M = NUM_ELITE * MaxThreads;
    for (i32 R = M-1; R > 0; R--)
    {
        i32 L = rand() % (R+1);
        memcpy(Buf, GetIndividuo(Elite, N, L), sizeof(i32)*N);
        memcpy(GetIndividuo(Elite, N, L), GetIndividuo(Elite, N, R), sizeof(i32)*N);
        memcpy(GetIndividuo(Elite, N, R), Buf, sizeof(i32)*N);
        r32 Tmpf = PuntuacionesElite[L];
        PuntuacionesElite[L] = PuntuacionesElite[R];
        PuntuacionesElite[R] = Tmpf;
    }
}

internal void
EjecutaSecuencial(i32 N,
        i32 TamPoblacion,
        rng *__restrict__ Rng,
        i32 *__restrict__ Poblacion,
        r32 *__restrict__ Puntuaciones,
        i32 *__restrict__ NuevaGeneracion,
        r32 *__restrict__ NuevasPuntuaciones,
        i32 *__restrict__ Elite,
        i32 *__restrict__ EliteNueva,
        r32 *__restrict__ PuntuacionesElite,
        i32 *__restrict__ out_Permutation,
        r32 *__restrict__ out_PermutationFitness,
        tsp_instance *__restrict__ Tsp)
{
    // Para cada elemento de la población elegimos 2 padres. Su posición de
    // generación anterior y el ganador de un torneo entre 5 elecciones aleatorias.
    // Se incluye elitismo, haciendo un torneo entre el ganador del torneo anterior y
    // una solución buena.
    for (i32 J = 0; J < TamPoblacion; J++)
    {
        i32 *Padre1 = GetIndividuo(Poblacion, N, J);
        i32 PosPadre = Torneo(Puntuaciones, TamPoblacion, Rng);
        i32 *Padre2 = GetIndividuo(Poblacion, N, PosPadre);
        if (Rng->Rng_r32() <= PROBA_ELITISMO)
        {
            i32 PosElite = Rng->Rng_i32(NUM_ELITE);
            i32 *PadreElite = GetIndividuo(Elite, N, PosElite);
            Padre2 = PadreElite;
        }
        i32 *Hijo = GetIndividuo(NuevaGeneracion, N, J);
        Recombina(Hijo, Padre1, Padre2, MUTACION, N, Rng);
    }
    // Evaluamos la población recién generada.
    for (i32 K = 0; K < TamPoblacion; K++)
    {
        i32 *Individuo = GetIndividuo(NuevaGeneracion, N, K);
        NuevasPuntuaciones[K] = Aptitud(N, Individuo, Tsp->Coords);
    }
    // Si el padre en la misma posición de un individuo era mejor, lo pasamos
    // a la siguiente generación con alguna probabilidad
    for (i32 K = 0; K < TamPoblacion; K++)
    {
        if (Puntuaciones[K] < NuevasPuntuaciones[K]
                && Rng->Rng_r32() > PROBA_HIJO_GRATIS)
        {
            i32 *Hijo = GetIndividuo(NuevaGeneracion, N, K);
            i32 *Padre = GetIndividuo(Poblacion, N, K);
            memcpy(Hijo, Padre, sizeof(i32)*N);
            NuevasPuntuaciones[K] = Puntuaciones[K];
        }
    }
    // Elegimos la mejor solución (la que minimiza la distancia) y la guardamos.
    i32 *Sol = GetIndividuo(NuevaGeneracion, N,
            ConsigueMejor(NuevasPuntuaciones, TamPoblacion));
    // Guardamos el mejor visto
    memcpy(out_Permutation, Sol, sizeof(i32)*N);
    // La nueva generación se vuelve la generación anterior para la siguiente iteración.
    std::swap(Poblacion, NuevaGeneracion);
    std::swap(Puntuaciones, NuevasPuntuaciones);

    ActualizaElite(Elite, EliteNueva, Poblacion, Puntuaciones, TamPoblacion, N, false);
    std::swap(Elite, EliteNueva);

    // Actualizamos las puntuaciones de la élite.
    for (i32 I = 0; I < NUM_ELITE; I++)
    {
        i32 *Individuo = GetIndividuo(Elite, N, I);
        PuntuacionesElite[I] = Aptitud(N, Individuo, Tsp->Coords);
    }
    // Guardamos el puntaje de la solución con mejor puntaje
    *out_PermutationFitness = Aptitud(N, out_Permutation, Tsp->Coords);
}

/**
 * Executes the euclidean TSP solver.
 *
 * @param Tsp A pointer to the TSP instance to be read from.
 * @param out_Permutation An array with N spaces, to be 0 through N-1.
 * @param Iterations Number or iteration to be executed.
 * 
 * @return 1 (everything ok).
 */
internal b32
Main(tsp_instance *__restrict__ Tsp,
        i32 *__restrict__ out_Permutation,
        u64 *__restrict__ Iterations,
        r32 Cutoff) 
{
    // Obtenemos los parámetros de la población
    i32 N = Tsp->N;
    i32 TamPoblacion = std::max(2*N, 1000*MaxThreads);

    size_t RngLen = sizeof(rng)*MaxThreads;
    size_t IslasLen = sizeof(i32) * N*TamPoblacion*(size_t)MaxThreads;
    size_t PuntuacionesTLen = sizeof(r32) * TamPoblacion*(size_t)MaxThreads;
    size_t IslasEliteLen = sizeof(i32) * NUM_ELITE*MaxThreads*(size_t)N;
    size_t PuntuacionEliteLen = sizeof(i32) * NUM_ELITE*(size_t)MaxThreads;
    size_t SolsPorHiloLen = sizeof(i32) * N*(size_t)MaxThreads;
    size_t FitPorHiloLen = sizeof(r32) * MaxThreads;

    size_t ArenaLen = RngLen + IslasLen +  PuntuacionesTLen +
        IslasLen + PuntuacionesTLen + IslasEliteLen + IslasEliteLen +
        PuntuacionEliteLen + SolsPorHiloLen + FitPorHiloLen;

    void *Arena = mmap(0, ArenaLen,
            PROT_READ | PROT_WRITE,
            MAP_PRIVATE | MAP_ANONYMOUS,
            -1, 0);

    if (Arena == MAP_FAILED)
    {
        IGNORE_RESULT(write(1, ERR_MALLOC, sizeof(ERR_MALLOC)));
        return 0;
    }

    rng *Rngs = (rng *)Arena;
    i32 *Islas = (i32*)(((u8*)Rngs) + RngLen);
    r32 *PuntuacionesT = (r32*)(((u8*)Islas) + IslasLen);
    i32 *NuevasIslas = (i32*)(((u8*)PuntuacionesT) + PuntuacionesTLen);
    r32 *NuevasPuntuacionesT = (r32*)(((u8*)NuevasIslas) + IslasLen);
    i32 *IslasElite = (i32*)(((u8*)NuevasPuntuacionesT) + PuntuacionesTLen);
    i32 *IslasEliteNuevas = (i32*)(((u8*)IslasElite) + IslasEliteLen);
    r32 *PuntuacionesEliteT = (r32*)(((u8*)IslasEliteNuevas) + IslasEliteLen);
    i32 *MejoresSols = (i32*)(((u8*)PuntuacionesEliteT) + PuntuacionEliteLen);
    r32 *MejoresCostos = (r32*)(((u8*)MejoresSols) + SolsPorHiloLen);

    #pragma omp parallel for
    for (i32 Idx = 0; Idx < MaxThreads; Idx++)
    {
        new (&Rngs[Idx]) rng();
        i32 *Poblacion = GetPoblacion(Islas, TamPoblacion, N, Idx);
        rng *Rng = Rngs + Idx;
        IniciaPoblacion(Poblacion, N, TamPoblacion, Rng);

        // Evalúa cada candidato
        r32 *Puntuaciones = GetDatoPoblacion(PuntuacionesT, TamPoblacion, Idx);
        for (i32 I = 0; I < TamPoblacion; I++)
        {
            i32 *Individuo = GetIndividuo(Poblacion, N, I);
            Puntuaciones[I] = Aptitud(N, Individuo, Tsp->Coords);
        }

        i32 *Elite = GetPoblacion(IslasElite, NUM_ELITE, N, Idx);
        i32 *EliteNueva = GetPoblacion(IslasEliteNuevas, NUM_ELITE, N, Idx);
        ActualizaElite(Elite, EliteNueva, Poblacion, Puntuaciones,
                TamPoblacion, N, true);

        r32 *PuntuacionesElite = GetDatoPoblacion(PuntuacionesEliteT,
                NUM_ELITE, Idx);
        for (i32 I = 0; I < NUM_ELITE; I++)
        {
            i32 *Individuo = GetIndividuo(Elite, N, I);
            PuntuacionesElite[I] = Aptitud(N, Individuo, Tsp->Coords);
        }
    }

    std::swap(IslasElite, IslasEliteNuevas);
    RemezclaElite(N, IslasElite, PuntuacionesEliteT);

    u64 MaxIt = *Iterations;
    // Ejecutamos el algoritmo el número de iteraciones especificadas.
    for (u64 I = 0; I < MaxIt; I++)
    {
        #pragma omp parallel for
        for (i32 Idx = 0; Idx < MaxThreads; Idx++)
        {
            rng *Rng = Rngs + Idx;
            i32 *Poblacion = GetPoblacion(Islas, TamPoblacion, N, Idx);
            r32 *Puntuaciones = GetDatoPoblacion(PuntuacionesT, TamPoblacion, Idx);
            i32 *NuevaGeneracion = GetPoblacion(NuevasIslas, TamPoblacion, N, Idx);
            r32 *NuevasPuntuaciones = GetDatoPoblacion(NuevasPuntuacionesT,
                    TamPoblacion, Idx);
            i32 *Elite = GetPoblacion(IslasElite, NUM_ELITE, N, Idx);
            i32 *EliteNuevas = GetPoblacion(IslasEliteNuevas, NUM_ELITE, N, Idx);
            r32 *PuntuacionesElite = GetDatoPoblacion(PuntuacionesEliteT, NUM_ELITE, Idx);
            i32 *out_Perm = MejoresSols + N*Idx;
            r32 *out_PermFit = MejoresCostos + Idx;
            EjecutaSecuencial(N, TamPoblacion, Rng, Poblacion, Puntuaciones,
                    NuevaGeneracion, NuevasPuntuaciones, Elite, EliteNuevas,
                    PuntuacionesElite, out_Perm, out_PermFit, Tsp);
        }

        i32 MejorS = ConsigueMejor(MejoresCostos, MaxThreads);
        i32 *S = GetIndividuo(MejoresSols, N, MejorS);
        memcpy(out_Permutation, S, sizeof(i32)*N);

        if (MejoresCostos[MejorS] <= Cutoff)
        {
            *Iterations = I+1;
            break;
        }

        std::swap(Islas, NuevasIslas);
        std::swap(PuntuacionesT, NuevasPuntuacionesT);
        std::swap(IslasElite, IslasEliteNuevas);
        RemezclaElite(N, IslasElite, PuntuacionesEliteT);
    }

    munmap(Arena, ArenaLen);

    return 1;
}

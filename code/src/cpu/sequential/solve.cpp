#pragma once

#include "defines.h"
#include "tsp.h"

#include <unistd.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <immintrin.h>
#include <algorithm>

#define ERR_MALLOC "No se pudo alojar suficiente memoria.\n"
#define NUM_ELITE 5
#define MUTACION 0.1f
#define PROBA_ELITISMO 0.1f
#define PROBA_HIJO_GRATIS 0.1f

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
 * Crea una población aleatoria de P soluciones de tamaño N.
 *
 * @param Poblacion El buffer con los N*P i32 a llenar.
 * @param N Tamaño de las soluciones o individuos.
 * @param P Número de individuos de la población.
 */
internal void
IniciaPoblacion(i32 *Poblacion, i32 N, i32 P)
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
            i32 Left = rand() % (Right + 1);
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
Torneo(r32 Puntuaciones[], i32 Len)
{
    i32 Random1 = rand() % Len;
    i32 Random2 = rand() % Len;
    i32 Random3 = rand() % Len;
    i32 Random4 = rand() % Len;
    i32 Random5 = rand() % Len;
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
Recombina(i32 *NuevoIndividuo, i32 *Padre1, i32 *Padre2, r32 Mutacion, i32 N)
{
    //Combinamos a los padres.
    i32* Solucion = NuevoIndividuo;
    std::vector<bool> Bitmap(N, 0);
    i32 Pos1 = rand() % N;
    i32 Pos2 = rand() % N;
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
    if (((r32)rand() / (RAND_MAX)) <= Mutacion)
    {
        i32 Random1 = rand() % N;
        i32 Random2 = rand() % N;
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
    // Agregamos semilla arbitraria al generador de números aleatorios
    srand(__rdtsc());

    // Creamos una población inicial aleatoria y la evaluamos.
    i32 N = Tsp->N;
    i32 TamPoblacion = std::max(2*N, 1000);

    i32 *Poblacion = (i32*)malloc(sizeof(i32) * N*TamPoblacion);
    
    // Si falla malloc, termina
    if (!Poblacion)
    {
        IGNORE_RESULT(write(1, ERR_MALLOC, sizeof(ERR_MALLOC)));
        return 0;
    }

    IniciaPoblacion(Poblacion, N, TamPoblacion);

    r32 *Puntuaciones = (r32*)malloc(sizeof(r32) * TamPoblacion);
    // Si falla malloc, termina
    if (!Puntuaciones)
    {
        free(Poblacion);
        IGNORE_RESULT(write(1, ERR_MALLOC, sizeof(ERR_MALLOC)));
        return 0;
    }

    // Evalúa cada candidato
    for (i32 I = 0; I < TamPoblacion; I++)
    {
        i32 *Individuo = GetIndividuo(Poblacion, N, I);
        Puntuaciones[I] = Aptitud(N, Individuo, Tsp->Coords);
    }

    // Creamos arreglos auxiliares para guardar la generación actual y la
    // anterior, o en este caso la generación anterior y la nueva.
    i32 *NuevaGeneracion = (i32*)malloc(sizeof(i32) * N*TamPoblacion);
    if (!NuevaGeneracion)
    {
        free(Poblacion);
        free(Puntuaciones);
        IGNORE_RESULT(write(1, ERR_MALLOC, sizeof(ERR_MALLOC)));
        return 0;
    }
    r32 *NuevasPuntuaciones = (r32*)malloc(sizeof(r32) * TamPoblacion);
    if (!NuevasPuntuaciones)
    {
        free(Poblacion);
        free(Puntuaciones);
        free(NuevaGeneracion);
        IGNORE_RESULT(write(1, ERR_MALLOC, sizeof(ERR_MALLOC)));
        return 0;
    }

    // Guardamos las (NUM_ELITE) mejores soluciones encontradas hasta el momento.
    i32 *Elite = (i32*)malloc(sizeof(i32) * NUM_ELITE*N);
    if (!Elite)
    {
        free(Poblacion);
        free(Puntuaciones);
        free(NuevaGeneracion);
        free(NuevasPuntuaciones);
        IGNORE_RESULT(write(1, ERR_MALLOC, sizeof(ERR_MALLOC)));
        return 0;
    }
    i32 *EliteNueva = (i32*)malloc(sizeof(i32) * NUM_ELITE*N);
    if (!EliteNueva)
    {
        free(Poblacion);
        free(Puntuaciones);
        free(NuevaGeneracion);
        free(NuevasPuntuaciones);
        free(Elite);
        IGNORE_RESULT(write(1, ERR_MALLOC, sizeof(ERR_MALLOC)));
        return 0;
    }
    r32 *PuntuacionesElite = (r32*)malloc(sizeof(r32) * NUM_ELITE);
    if (!PuntuacionesElite)
    {
        free(Poblacion);
        free(Puntuaciones);
        free(NuevaGeneracion);
        free(NuevasPuntuaciones);
        free(Elite);
        free(EliteNueva);
        IGNORE_RESULT(write(1, ERR_MALLOC, sizeof(ERR_MALLOC)));
        return 0;
    }
    ActualizaElite(Elite, EliteNueva, Poblacion, Puntuaciones, TamPoblacion, N, true);
    std::swap(Elite, EliteNueva);

    for (i32 I = 0; I < NUM_ELITE; I++)
    {
        i32 *Individuo = GetIndividuo(Elite, N, I);
        PuntuacionesElite[I] = Aptitud(N, Individuo, Tsp->Coords);
    }

    u64 MaxIt = *Iterations;
    // Ejecutamos el algoritmo el número de iteraciones especificadas.
    for (u64 I = 0; I < MaxIt; I++)
    {
        // Para cada elemento de la población elegimos 2 padres. Su posición de
        // generación anterior y el ganador de un torneo entre 5 elecciones aleatorias.
        // Se incluye elitismo, haciendo un torneo entre el ganador del torneo anterior y
        // una solución buena.
        for (i32 J = 0; J < TamPoblacion; J++)
        {
            i32 *Padre1 = GetIndividuo(Poblacion, N, J);
            i32 PosPadre = Torneo(Puntuaciones, TamPoblacion);
            i32 *Padre2 = GetIndividuo(Poblacion, N, PosPadre);
            if (((r32)rand() / (RAND_MAX)) <= PROBA_ELITISMO)
            {
                i32 PosElite = rand() % NUM_ELITE;
                i32 *PadreElite = GetIndividuo(Elite, N, PosElite);
                Padre2 = PadreElite;
            }
            i32 *Hijo = GetIndividuo(NuevaGeneracion, N, J);
            Recombina(Hijo, Padre1, Padre2, MUTACION, N);
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
                    && ((r32)rand() / (RAND_MAX)) > PROBA_HIJO_GRATIS)
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
        // Revisamos que la solución con mejor puntaje al momento sea suficiente.
        if(Aptitud(N, out_Permutation, Tsp->Coords) <= Cutoff)
        {
            *Iterations = I + 1;
            break;
        }
    }

    free(Poblacion);
    free(Puntuaciones);
    free(NuevaGeneracion);
    free(NuevasPuntuaciones);
    free(Elite);
    free(EliteNueva);
    free(PuntuacionesElite);

    return 1;
}

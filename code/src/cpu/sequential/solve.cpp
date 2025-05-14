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

/**
 * Crea una población aleatoria de P soluciones de tamaño N.
 *
 * @param N Tamaño de las soluciones o individuos.
 * @param P Número de individuos de la población.
 * 
 * @return Poblacion - La población aleatoria formada.
 */
i32**
CrearPoblacion(i32 N, i32 P)
{
    i32 **Poblacion = (i32**)malloc(sizeof(i32*) * P);
    for (i32 I = 0; I < P; I++)
    {
        i32 *Arr = (i32*)malloc(sizeof(i32) * N);
        // Llenamos el arreglo con los números del 0 a N-1.
        for (i32 J = 0; J < N; J++)
        {
            Arr[J] = J;
        }
        // Hacemos shuffle de los números.
        for (i32 Right = N-1; Right > 0; Right--)
        {
            i32 Left = rand() % (Right + 1);
            i32 Tmp = Arr[Left];
            Arr[Left] = Arr[Right];
            Arr[Right] = Tmp;
        }
        Poblacion[I] = Arr;
    }
    return Poblacion;
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
r32
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
r32
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
i32
ConsigueMejor(r32 Puntuaciones[], i32 Len)
{
    r32 menor = Puntuaciones[0];
    i32 posicion = 0;
    for (i32 I = 1; I < Len; I++)
    {
        if (Puntuaciones[I] < menor)
        {
            menor = Puntuaciones[I];
            posicion = I;
        }
    }
    return posicion;
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
i32
Torneo(r32 Puntuaciones[], i32 Len)
{
    i32 Random1 = rand() % Len;
    i32 Random2 = rand() % Len;
    i32 Random3 = rand() % Len;
    i32 Random4 = rand() % Len;
    i32 Pos[] = {Random1, Random2, Random3, Random4};
    r32 Arr[] = {Puntuaciones[Random1], Puntuaciones[Random2], Puntuaciones[Random3], Puntuaciones[Random4]};
    return Pos[ConsigueMejor(Arr, 5)];
}

/**
 * Toma dos soluciones padre y las combina para generar una nueva solución hija.
 *
 * @param Padre1 Solución del primer padre, la cual se recombinará.
 * @param Padre2 Solución del segundo padre, la cual se recombinará.
 * @param Mutacion Probabilidad de que la solución mute.
 * @param N Número de ciudades.
 * 
 * @return Solucion - Solución que combina a sus padres. Puede tener mutaciones.
 */
i32*
Recombina(i32 *Padre1, i32 *Padre2, r32 Mutacion, i32 N)
{
    //Combinamos a los padres.
    i32* Solucion = (i32*)malloc(sizeof(i32) * N);
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
    
    // Aplicamos mutación itercambiando dos ciudades aleatorias.
    if (((r32)rand() / (RAND_MAX)) <= Mutacion)
    {
        i32 Random1 = rand() % N;
        i32 Temp = Solucion[Random1];
        i32 Random2 = rand() % N;
        Solucion[Random1] = Solucion[Random2];
        Solucion[Random2] = Temp;
    }
    return Solucion;
}

/**
 * Encuentra los mejores 5 valores entre la generación y las respuestas
 * previamente contenidas en la élite.
 *
 * @param Elite Un arreglo que contiene las 5 mejores soluciones encontradas al
 *              momento.
 * @param Poblacion Arreglo con las soluciones de la generación actual.
 * @param Puntuaciones Puntos asignados a cada solución de la generación actual.
 * 
 * @return Nuevo arreglo con los 5 mejores individuos encontrados.
 */
i32**
ActualizaElite(i32 **Elite, i32 **Poblacion, r32 *Puntuaciones, i32 P, i32 N)
{
    // Guardamos la puntuación y el índice de la solución en Poblacion.
    std::vector<std::pair<r32, i32>> Pares;
    for (i32 I = 0; I < P; ++I) {
        Pares.push_back({Puntuaciones[I], I});
    }

    // Ordenamos las puntuaciones en órden ascendente
    //(dando prioridad a las menores puntuaciones).
    std::sort(Pares.begin(), Pares.end());

    // Guardamos todos los individuos, tanto de la élite como de la población.
    std::vector<i32*> CandidatosAll;
    std::vector<r32> PuntuacionesAll;
    
    for (i32 I = 0; I < 5; ++I) {
        if (Elite[I] != nullptr) {
            CandidatosAll.push_back(Elite[I]);
            PuntuacionesAll.push_back(std::numeric_limits<r32>::max());
        }
    }

    i32 Candidatos = std::min(P, 5);
    for (i32 I = 0; I < Candidatos; ++I) {
        CandidatosAll.push_back(Poblacion[Pares[I].second]);
        PuntuacionesAll.push_back(Pares[I].first);
    }

    // Creamos un vector de pares para ordenar todos los candidatos.
    std::vector<std::pair<r32, i32*>> CandidatosOrdenados;
    for (size_t I = 0; I < CandidatosAll.size(); ++I) {
        CandidatosOrdenados.push_back({PuntuacionesAll[I], CandidatosAll[I]});
    }

    // Ordenamos a los candidatos.
    std::sort(CandidatosOrdenados.begin(), CandidatosOrdenados.end());

    // Crear el nuevo arreglo de la élite
    i32** EliteNueva = (i32**)malloc(sizeof(i32*) * 5);

    // Copiamos las 5 mejores soluciones a la nueva élite.
    for (i32 I = 0; I < 5; ++I) {
        EliteNueva[I] = (i32*)malloc(sizeof(i32) * N);
        for (i32 J = 0; J < N; ++J) {
            EliteNueva[I][J] = CandidatosOrdenados[I].second[J];
        }
    }
    return EliteNueva;
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
b32
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
    i32 **Poblacion = CrearPoblacion(N, TamPoblacion);
    r32 *Puntuaciones = (r32*)malloc(sizeof(r32) * TamPoblacion);
    for (i32 I = 0; I < TamPoblacion; I++)
    {
        Puntuaciones[I] = Aptitud(N, Poblacion[I], Tsp->Coords);
    }

    // Creamos arreglos auxiliares para guardar la generación actual y la anterior.
    // o en este caso la generación anterior y la nueva.
    i32 **NuevaGeneracion = (i32**)malloc(sizeof(i32*) * TamPoblacion);
    r32 *NuevasPuntuaciones = (r32*)malloc(sizeof(r32) * TamPoblacion);

    // Guardamos las 5 mejores soluciones encontradas hasta el momento.
    i32 **Elite = (i32**)malloc(sizeof(i32*) * 5);
    r32 *PuntuacionesElite = (r32*)malloc(sizeof(r32) * 5);

    Elite = ActualizaElite(Elite, Poblacion, Puntuaciones, TamPoblacion, N);

    for (i32 I = 0; I < 5; I++)
    {
        PuntuacionesElite[I] = Aptitud(N, Elite[I], Tsp->Coords);
    }

    // Ejecutamos el algoritmo el número de iteraciones especificadas.
    for (u64 I = 0; I < *Iterations; I++)
    {
        // Para cada elemento de la población elegimos 2 padres. Su posición de
        // generación anterior y el ganador de un torneo entre 5 elecciones aleatorias.
        // Se incluye elitismo, haciendo un torneo entre el ganador del torneo anterior y
        // una solución buena.
        for (i32 J = 0; J < TamPoblacion; J++)
        {
            i32 *Padre1 = Poblacion[J];
            i32 PosElite = rand() % 5;
            i32 *PadreElite = Elite[PosElite];
            i32 PosPadre = Torneo(Puntuaciones, TamPoblacion);
            i32 *Padre2 = Poblacion[PosPadre];
            if (PuntuacionesElite[PosElite] < Puntuaciones[PosPadre])
            {
                Padre2 = PadreElite;
            }
            NuevaGeneracion[J] = Recombina(Padre1, Padre2, 0.5f, N);
        }
        // Evaluamos la población recién generada.
        for (i32 K = 0; K < TamPoblacion; K++)
        {
            NuevasPuntuaciones[K] = Aptitud(N, NuevaGeneracion[K], Tsp->Coords);
        }
        // Elegimos la mejor solución (la que minimiza la distancia) y la guardamos.
        i32 *Sol = NuevaGeneracion[ConsigueMejor(NuevasPuntuaciones, TamPoblacion)];
        memcpy(out_Permutation,Sol,sizeof(i32)*N);
        // La nueva generación se vuelve la generación anterior para la siguiente iteración.
        memcpy(Poblacion,NuevaGeneracion,sizeof(i32*)*TamPoblacion);
        memcpy(Puntuaciones,NuevasPuntuaciones,sizeof(r32)*TamPoblacion);

        Elite = ActualizaElite(Elite, Poblacion, Puntuaciones, TamPoblacion, N);

        // Actualizamos las puntuaciones de la élite.
        for (i32 I = 0; I < 5; I++)
        {
            PuntuacionesElite[I] = Aptitud(N, Elite[I], Tsp->Coords);
        }
        // Revisamos que la solución con mejor puntaje al momento sea suficiente.
        if(Aptitud(N, out_Permutation, Tsp->Coords) <= Cutoff)
        {
            *Iterations = I;
            break;
        }
    }
    return 1;
}
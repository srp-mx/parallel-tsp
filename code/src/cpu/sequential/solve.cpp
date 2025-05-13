#pragma once

#include "defines.h"
#include "tsp.h"

#include <unistd.h>
#include <stdio.h>
#include <malloc.h>
#include <bits/stdc++.h>
#include <vector>

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
    i32 **Poblacion = (i32**)malloc(sizeof(int*) * P);
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
    i32 menor = Puntuaciones[0];
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
    i32 Random5 = rand() % Len;
    i32 Pos[] = {Random1, Random2, Random3, Random4, Random5};
    r32 Arr[] = {Puntuaciones[Random1], Puntuaciones[Random2], Puntuaciones[Random3], Puntuaciones[Random4], Puntuaciones[Random5]};
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
    // Creamos una población inicial aleatoria y la evaluamos.
    i32 N = Tsp->N;
    i32 **Poblacion = CrearPoblacion(N, 2*N);
    r32 *Puntuaciones = (r32*)malloc(sizeof(r32) * N);
    for (i32 I = 0; I < N; I++)
    {
        Puntuaciones[I] = Aptitud(N, Poblacion[I], Tsp->Coords);
    }

    // Creamos arreglos auxiliares para guardar la generación actual y la anterior.
    // o en este caso la generación anterior y la nueva.
    i32 **NuevaGeneracion = (i32**)malloc(sizeof(i32)*2*N*N);
    r32 *NuevasPuntuaciones = (r32*)malloc(sizeof(r32) * N);

    // Ejecutamos el algoritmo el número de iteraciones especificadas.
    for (u64 I = 0; I < *Iterations; I++)
    {
        // Para cada elemento de la población elegimos 2 padres. Su posición de
        // generación anterior y el ganador de un torneo entre 5 elecciones aleatorias.
        for (i32 J = 0; J < N; J++)
        {
            i32 *Padre1 = Poblacion[J];
            i32 *Padre2 = Poblacion[Torneo(Puntuaciones, 2*N)];
            NuevaGeneracion[J] = Recombina(Padre1, Padre2, 0.5f, N);
        }
        // Evaluamos la población recién generada.
        for (i32 K = 0; K < N; K++)
        {
            NuevasPuntuaciones[K] = Aptitud(N, Poblacion[K], Tsp->Coords);
        }
        // Elegimos la mejor solución (la que minimiza la distancia) y la guardamos.
        out_Permutation = NuevaGeneracion[ConsigueMejor(NuevasPuntuaciones, N)];
        // La nueva generación se vuelve la generación anterior para la siguiente iteración.
        Poblacion = NuevaGeneracion;
        Puntuaciones = NuevasPuntuaciones;
        // Revisamos que la solución con mejor puntaje al momento sea suficiente.
        if(Aptitud(N, out_Permutation, Tsp->Coords) <= Cutoff)
        {
            *Iterations = I;
            break;
        }
    }
    return 1;
}
#pragma once

#include "defines.h"
#include "tsp.h"

#include <unistd.h>
#include <stdio.h>
#include <malloc.h>
#include <bits/stdc++.h>

/**
 * 1.Crear una población inicial de tamaño n.
 * 2.Calcular el valor de aptitud para cada individuo.
 * 3.Seleccionar un subconjunto de soluciones candidatas para formar parte del
 * conjunto de apareamiento.
 * 4.Realizar la recombinación entre los cromosomas padres.
 * 5.Aplicar mutación.
 * 6.Repetir hasta generar el número de generaciones establecidas
 */

i32**
CrearPoblacion(i32 N, i32 P)
{
    i32 *Arr = (i32*)malloc(sizeof(i32) * N);
    i32 **Poblacion = (i32**)malloc(sizeof(Arr) * P);
    for (i32 I = 0; I < P; I++)
    {
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

i32
Aptitud(i32 *Individuo, v2 Coords[])
{

}

void
SeleccionaPadres() {}

void
Recombina() {}

void
Mutacion() {}

b32
Main(tsp_instance *__restrict__ Tsp,
        i32 *__restrict__ out_Permutation,
        u64 *__restrict__ Iterations,
        r32 Cutoff) 
{
    i32 N = Tsp->N;
    i32 **Poblacion = CrearPoblacion(N, 2*N);
    i32 *Puntuaciones = (i32*)malloc(sizeof(i32) * N);
    for (i32 I = 0; I < N; I++)
    {
        Puntuaciones[I] = Aptitud(Poblacion[I], Tsp->Coords);
    }
    return 1;
}

/**
 * Tsp->Coords[u]
0 <= u < Tsp->N
dx = coordA.X - coordB.X
dy = coordA.Y - coordB.Y
dist = sqrtf(dx*dx + dy*dy)

acc = 0.0f;
for(i = 0; i < n; i++)
{
acc += dist(tu_perm[i], tu_perm[(i+1)%n])
}


100 200 300
1/6 2/6 3/6
5/6 4/6 3/6
100 200 300


(600-100)/600
1/6 2/6 3/6
0.34


0.34-1/6 mas o menos 0.16


Poblacion[i]
Poblacion[j]
u32 Poblacion[P][N] = {};


u32 PoblacionNueva[P][N] = {};
for (i = 0; i < N; i++)
{

}
Padre1 = Poblacion[i]


Padre2 = Poblacion[torneo()]


code/data/main
solver cpuseq
b Run
[tsp]$ run
solve_Solve()


for (i32 Right = N-1; Right > 0; Right--)
    {
        i32 Left = curand(&Rng) % (Right + 1);
        i32 Tmp = Entry[Left];
        Entry[Left] = Entry[Right];
        Entry[Right] = Tmp;
    }
 */
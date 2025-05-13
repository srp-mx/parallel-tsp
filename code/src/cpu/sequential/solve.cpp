#pragma once

#include "defines.h"
#include "tsp.h"

#include <unistd.h>
#include <stdio.h>
#include <malloc.h>
#include <bits/stdc++.h>

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

r32
Distancia(i32 Ciudad1, i32 Ciudad2, v2 Coords[])
{
    r32 Dx = Coords[Ciudad1].X - Coords[Ciudad2].X;
    r32 Dy = Coords[Ciudad1].Y - Coords[Ciudad2].Y;
    return sqrtf(Dx*Dx + Dy*Dy);
}

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

i32*
Recombina(i32 *Padre1, i32 *Padre2, r32 Mutacion, i32 N)
{
    //Combinamos a los padres.
    // TODO: combinacion de soluciones
    i32* Solucion = 0;
    // Aplicamos mutación.
    if (((r32)rand() / (RAND_MAX)) <= Mutacion)
    {
        i32 Random1 = rand() % N;
        i32 Temp = Solucion[Random1];
        i32 Random2 = rand() % N;
        Solucion[Random1] = Solucion[Random2];
        Solucion[Random2] = Temp;
    }
}

b32
Main(tsp_instance *__restrict__ Tsp,
        i32 *__restrict__ out_Permutation,
        u64 *__restrict__ Iterations,
        r32 Cutoff) 
{
    i32 N = Tsp->N;
    i32 **Poblacion = CrearPoblacion(N, 2*N);
    r32 *Puntuaciones = (r32*)malloc(sizeof(r32) * N);
    for (i32 I = 0; I < N; I++)
    {
        Puntuaciones[I] = Aptitud(N, Poblacion[I], Tsp->Coords);
    }

    i32 **NuevaGeneracion = (i32**)malloc(sizeof(i32)*2*N*N);
    r32 *NuevasPuntuaciones = (r32*)malloc(sizeof(r32) * N);

    for (u64 I = 0; I < *Iterations; I++)
    {
        for (i32 J = 0; J < N; J++)
        {
            i32 *Padre1 = Poblacion[J];
            i32 *Padre2 = Poblacion[Torneo(Puntuaciones, 2*N)];
            NuevaGeneracion[J] = Recombina(Padre1, Padre2, 0.5f, N);
        }
        out_Permutation = NuevaGeneracion[ConsigueMejor(NuevasPuntuaciones, N)];
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

/**
 * 
100 200 300
1/6 2/6 3/6
5/6 4/6 3/6
100 200 300


(600-100)/600
1/6 2/6 3/6
0.34


0.34-1/6 mas o menos 0.16


code/data/main
solver cpuseq
b Run
[tsp]$ run
solve_Solve()




bool repetido = false;
i32 random = rand() % N;
for (i32 k = 0; k < j; k++)
{
    if (random == arr[k])
    {
        repetido = true;
        break;
    }
}
if (!repetido)
{
    arr[j] = random;
    j++;
}
 */
#pragma once

extern "C"
{
    /**
     * Vector 2D.
     */
    typedef struct v2
    {
        r32 X, Y;
    } v2;

    /**
     * TSP instance data.
     */
    typedef struct tsp_instance
    {
        i32 N;
        v2 Coords[];
    } tsp_instance;
}

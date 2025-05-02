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

#include <stdlib.h>

#include "defines.h"
#include "tsp.h"

/**
 * Auxiliary data to parse the TSP instance file.
 */
struct parsed_tsp
{
    char *Name, *Comment;
    tsp_instance *Tsp;
    i32 N;
    b32 Ok;
};

/**
 * Moves the input string into the next line if it exists.
 *
 * @param Input A pointer to a char* with the text to parse.
 *
 * @return 1 if the input ended, 0 if it continues.
 */
internal b32
NextLine(char **Input)
{
    while (**Input && **Input != '\n')
    {
        (*Input)++;
    }
    if (!**Input)
    {
        return 1;
    }
    (*Input)++;
    if (**Input)
    {
        if ((*Input)[0] == 'E' && (*Input)[1] == 'O' && (*Input)[2] == 'F')
        {
            return 0;
        }
        return 0;
    }
    return 1;
}

/**
 * Advances a string into the first whitespace, null or ':' position after this
 * one, makes it null, and finally advances the string one place if it was not
 * null before. That way we create a string from some unspecified location
 * (ideally the starting position of the string) up to where the "word" "ends".
 *
 * @param End A pointer to the string to modify.
 *
 * @return 0 if the input ended, -1 if it reached a newline, 1 otherwise.
 */
internal b32
Tokenize(char **End)
{
    while (**End && **End != ' ' && **End != '\n' && **End != ':')
    {
        (*End)++;
    }
    if (**End)
    {
        **End = 0;
        (*End)++;
        if (**End == '\n')
        {
            return -1;
        }
        return 1;
    }
    return 0;
}

/**
 * Skips whitespace and newlines and parses what it reaches as an r32.
 *
 * @param Input A pointer to the input string, to be advanced.
 * @param out_Read A pointer to the r32 value to be assigned.
 *
 * @return If the parsing of the next number was successful.
 */
internal b32
ReadNextNum(char **Input, r32 *out_Read)
{
    while (**Input == ' ' || **Input == '\n')
    {
        (*Input)++;
    }
    if (!**Input) return 0;
    *out_Read = strtof(*Input, Input);
    return **Input == ' ' || **Input == '\n' || **Input == 0;
}

/**
 * Parses all the number triplets in the NODE_COORD_SECTION.
 *
 * @param N The number of coordinates.
 * @param Coords A pointer to the coordinate buffer.
 * @param Pointer to the first character in the section.
 *
 * @return If everything went ok.
 */
internal b32
ParseCoords(i32 N, v2 *Coords, char *Input)
{
    r32 Idx;
    for (i32 I = 0; I < N; I++, Coords++)
    {
        if (!ReadNextNum(&Input, &Idx)) return 0;
        if (!ReadNextNum(&Input, &(Coords->X))) return 0;
        if (!ReadNextNum(&Input, &(Coords->Y))) return 0;
    }
    return 1;
}

/**
 * Reads and skips or stores each label in the TSP file, such as NAME, COMMENT,
 * DIMENSION, NODE_COORD_SECTION, and so on.
 *
 * @param Parsed A pointer to the parsing data to be read and modified.
 * @param Behind Pointer to the string that stays behind to begin substrings.
 * @param Ahead Pointer to the string that stays ahead to add null terminators.
 *
 * @return If the reading must stop.
 */
internal b32
ReadDirective(parsed_tsp *Parsed, char **Behind, char **Ahead)
{
    if (!Tokenize(Ahead))
    {
        return 1;
    }

    while (**Ahead && (**Ahead == ' ' || **Ahead == ':'))
    {
        (*Ahead)++;
    }

    if (!strcmp(*Behind, "NAME"))
    {
        *Behind = *Ahead;
        if (!Tokenize(Ahead))
        {
            return 1;
        }
        Parsed->Name = *Behind;
    }
    else if (!strcmp(*Behind, "COMMENT"))
    {
        if (!Parsed->Comment)
        {
            Parsed->Comment = *Ahead;
        }
        while (**Ahead != '\n')
        {
            (*Ahead)++;
        }
        *((*Ahead)++) = 0;
    }
    else if (!strcmp(*Behind, "DIMENSION"))
    {
        Parsed->N = strtol(*Ahead, Ahead, 10);
        (*Ahead)++;
    }
    else if (!strcmp(*Behind, "NODE_COORD_SECTION"))
    {
        if (Parsed->N <= 0)
        {
            return 1;
        }
        Parsed->Tsp = (tsp_instance*)mmap(0,
                sizeof(tsp_instance) + sizeof(v2)*Parsed->N,
                PROT_READ | PROT_WRITE,
                MAP_PRIVATE | MAP_ANONYMOUS,
                -1, 0);
        Parsed->Tsp->N = Parsed->N;
        Parsed->Ok = ParseCoords(Parsed->N, Parsed->Tsp->Coords, *Ahead);
        return 1;
    }
    else
    {
        b32 Ends = NextLine(Ahead);
        *Behind = *Ahead;
        return Ends;
    }
    
    *Behind = *Ahead;
    return 0;
}

/**
 * Parses a string into a TSP parsing result.
 *
 * @param Input The text to parse.
 *
 * @return The TSP parsing result, which includes some metadata and the TSP
 *         instance.
 */
internal parsed_tsp
ParseTsp(char *Input)
{
    parsed_tsp Parsed = {};

    char *Behind, *Ahead;
    Behind = Ahead = Input;

    while (!ReadDirective(&Parsed, &Behind, &Ahead));

    return Parsed;
}

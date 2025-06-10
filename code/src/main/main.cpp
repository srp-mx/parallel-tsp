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

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <dlfcn.h>

#define MAIN_MAJOR_VERSION 2

#include "defines.h"
#include "config.h"
#include "utils.h"

#define INSTANT_WRITE(StringLiteral) write(1, StringLiteral, sizeof(StringLiteral))

#include "commands.h"

/**
 * Reads a command from the standard input.
 *
 * @return True if the parsing may continue, false if it ended.
 */
internal b32
ReadCommand(config *Config)
{
    char LineBuffer[MAX_FILENAME];
    char *FgetsStatus = fgets(LineBuffer, MAX_FILENAME-32, stdin);
    LineBuffer[MAX_FILENAME-32] = 0;
    if (feof(stdin))
    {
        INSTANT_WRITE("\nAdiós :)\n");
        return 0;
    }
    if (!FgetsStatus)
    {
        INSTANT_WRITE("Tu entrada es demasiado larga\n");
        return 1;
    }
    char *CommandStart, *CommandEnd, *Arg, *ArgEnd;
    for (CommandStart = LineBuffer; *CommandStart == ' '; CommandStart++);
    for (CommandEnd = CommandStart; *CommandEnd != ' ' && *CommandEnd != '\n' && *CommandEnd; CommandEnd++);
    for (Arg = CommandEnd; *Arg == ' '; Arg++);
    for (ArgEnd = Arg; *ArgEnd != '\n' && *ArgEnd; ArgEnd++);
    CommandEnd--;
    *ArgEnd = 0;
    command *Command = 0;
    for (size_t I = 0; I < CMD_COUNT; I++)
    {
        size_t CmpLen =
            Max(size_t(CommandEnd-CommandStart+1), COMMAND_NAME_LENS[I]);
        if (strncmp(CommandStart, COMMAND_NAMES[I], CmpLen))
            continue;
        Command = COMMANDS[I]; break;
    }
    if (Command) Command(Config, Arg);
    else INSTANT_WRITE(
            "No se reconoció el comando. "
            "Escribe help si tienes dudas.\n");
    return ContinueRepl;
}

/** main */
int
main()
{
    config Config = {};
    Config.Iterations = 1;
    Config.Executions = 1;
    do { INSTANT_WRITE("\033[1m[tsp]\033[0m$ "); } while (ReadCommand(&Config));

    if (Config.Solver.Handle)
    {
        dlclose(Config.Solver.Handle);
    }

    return 0;
}

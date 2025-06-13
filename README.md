# Parallel TSP
> This is a tool to use and benchmark different 2D Euclidean TSP solvers.

## Documentation

### Building

You will need `lualatex` and many LaTeX packages. If you're on Linux, we
recommend searching your distro's software repository for a package called
something akin to `texlive-all` or `texlive-full`.

Once inside each documentation's `latex` directory, if you're on Linux, you can
simply use
```bash
./run
```

### Index
- [Project proposal](docs/propuesta/latex)
- [Project report](docs/reporte/latex)

## Compilation

We assume you're on Linux, since we're not providing support for other
platforms at the moment.

You will need, at the very least, an `x86_64` CPU and the `g++` compiler,
oftentimes included in a package on your distro called `build-essential` or
`base-devel`.

### Entrypoint

You will first need to compile the program that allows you to setup the
execution configuration. This program is found under
[`code/src/main/`](code/src/main/).

You may modify the `build` script as you wish, and then execute the compilation
by running
```bash
./build
```

If you already have some solvers built, you may then execute the program using
```bash
./run
```

### Solvers

These are dynamic libraries that find the desired candidate solutions for the
2D Euclidean TSP instances. They are separate programs from the `main`
executable, since they may require special hardware and/or software support,
as well as being near-realtime debuggable by hot-swapping them on runtime.

That is, you can run the `./main` executable, and then recompile any solver
whenever you want to make a change without needing to rerun the program.
The only requirement is that you unload the library before doing so, by
swapping to another solver before compiling.

#### Dummy

The `dummy` solver simply returns the nodes in the order they are received
and is only for debugging purposes.

All you need to build it is to go into [`code/src/dummy/`](code/src/dummy/) and
run
```bash
./build
```

#### Sequential genetic

The `cpuseq` solver implements a genetic algorithm in sequential execution, and
it serves as a point of comparison.

All you need to build it is to go into
[`code/src/cpu/sequential/`](code/src/cpu/sequential/) and run
```bash
./build
```

#### Parallel genetic

The `cpupar` solver simply applies the fork-join technique into the same
algorithm as `cpuseq`. This serves as a point of comparison to a naïve parallel
implementation in CPU.

You will need to install `OpenMP` in your machine to compile it, but if you
already have the shared object file, you don't need it to run it.

To build it, go into
[`code/src/cpu/parallel/`](code/src/cpu/parallel/) and run
```bash
./build
```

#### GPU genetic

The `gpu` solver implements a genetic algorithm bespoke to the GPU's
architecture, and as such, requires special hardware and software.

To build it, you need to install the `CUDA Toolkit`, but you don't necessarily
need a GPU just yet. You only need to have the `nvcc` compiler available and
functional.

To execute it, you need a Pascal NVIDIA GPU or later, but if you already have
the shared object file, you don't need the `CUDA Toolkit` to run it.

To build it, go into
[`code/src/gpu/`](code/src/gpu/) and run
```bash
./build
```

## Execution

By default, once you build the entrypoint program and some solvers, their
respective ELF files will be found on [`code/data/`](code/data/). You can
execute `./main` from there, or as mentioned before, go into the
`code/src/main/` directory and execute `./run`.

Once inside the program, you can type the following for instructions
```
help
```

There are 10 commands, all of which admit either 0 or 1 argument:
- `help` Shows a help dialogue.
- `problem` Loads a problem in `code/data/tsp/*.tsp`.
- `solver` Loads a solver, unloading the one active beforehand.
- `iterations` Sets the maximum number of iterations per execution.
- `executions` Sets the number of executions from scratch.
- `cutoff` Sets the solution cost cutoff.
- `parallelism` Sets the amount of parallel "units" to use.
- `config` Shows the current configuration to launch the solver with.
- `exit` Exits the program, also doable with `ctrl-D`.
- `run` Executes the solver with the current configuration.

For example, `problem a280` will load the problem on the file
`code/data/tsp/a280.tsp`, while `solver dummy` will load the library
in `code/data/dummy.so`.

The difference between `iterations` and `executions` is that an iteration
generally depends on the last one, while each execution is a rerun of the
solver from scratch, completely oblivious to the last one. Doing multiple
executions allow us to get a mean and sample standard deviation, as well as
other statistical measurements, over multiple runs that share the same
configuration.

Once run, the program will display some execution information measured by
running the solvers, which will be saved in
[`code/data/logs/`](code/data/logs/) as a timestamped JSON file.

## Adding solvers

Build a dynamic library that includes, or at least conforms to, the file found
in [`code/src/include/solver.h`](code/src/include/solver.h). Make sure to
implement all functions.

Add your `.so` file, desirably with a short but unique name, to
[`code/data/`](code/data/).

That's it. It's as simple as that.

You can now load, for example, `my_solver.so` by simply using
```
solver my_solver
```

#!/usr/bin/env python3

from math import log10
import sys
import subprocess
from datetime import datetime

# Experiment data

selected_problems = [
    "berlin52",
    "kroB100",
    "lin105",
    "tsp225",
    "pcb442",
    "d657",
    "rat783",
    "pcb1173",
]
cutoff_percent = lambda prob_size : 1 + (0.01 * round(1+log10(prob_size)**2))
executions = 50
max_iterations = 250_000

# Information on some problems

# Maps problem to their known/reported lower bound.
# In cases where the lower bound was a range, we have taken the upper value.
lower_bounds = {
    "eil51":    426.0,
    "berlin52": 7542.0,
    "eil76":    538.0,
    "pr76":     108159.0,
    "rat99":    1211.0,
    "rd100":    7910.0,
    "kroA100":  21282.0,
    "kroB100":  22141.0,
    "kroC100":  20749.0,
    "kroD100":  21294.0,
    "kroE100":  22068.0,
    "eil101":   629.0,
    "lin105":   14379.0,
    "pr107":    44303.0,
    "pr124":    59030.0,
    "bier127":  118282.0,
    "pr136":    96772.0,
    "pr144":    58537.0,
    "ch150":    6528.0,
    "kroA150":  26524.0,
    "kroB150":  26130.0,
    "pr152":    73682.0,
    "u159":     42080.0,
    "rat195":   2323.0,
    "d198":     15780.0,
    "kroA200":  29368.0,
    "kroB200":  29437.0,
    "ts225":    126643.0,
    "tsp225":   3919.0,
    "pr226":    80369.0,
    "gil262":   2378.0,
    "pr264":    49135.0,
    "a280":     2579.0,
    "pr299":    48191.0,
    "lin318":   42029.0,
    "linhp318": 41345.0,
    "rd400":    15281.0,
    "fl417":    11861.0,
    "pr439":    107217.0,
    "pcb442":   50778.0,
    "d493":     35002.0,
    "u574":     36905.0,
    "rat575":   6773.0,
    "p654":     34643.0,
    "d657":     48912.0,
    "u724":     41910.0,
    "rat783":   8806.0,
    "pr1002":   259045.0,
    "u1060":    224094.0,
    "vm1084":   239297.0,
    "pcb1173":  56892.0,
    "d1291":    50801.0,
    "rl1304":   252948.0,
    "rl1323":   270199.0,
    "nrw1379":  56638.0,
    "fl1400":   20127.0,
    "u1432":    152970.0,
    "fl1577":   22249.0,
    "d1655":    62128.0,
    "vm1748":   336556.0,
    "u1817":    57201.0,
    "rl1889":   316536.0,
    "d2103":    80450.0,
    "u2152":    64253.0,
    "u2319":    234256.0,
    "pr2392":   378032.0,
    "pcb3038":  137694.0,
    "fl3795":   28772.0,
    "fnl4461":  182566.0,
    "rl5915":   565530.0,
    "rl5934":   556045.0,
    "rl11849":  923368.0,
    "usa13509": 19982889.0,
    "brd14051": 469445.0,
    "d15112":   1573152.0,
    "d18512":   645488.0
}

# Maps problem name with problem size
dimensions = {
    "eil51":    51,
    "berlin52": 52,
    "eil76":    76,
    "pr76":     76,
    "rat99":    99,
    "rd100":    100,
    "kroA100":  100,
    "kroB100":  100,
    "kroC100":  100,
    "kroD100":  100,
    "kroE100":  100,
    "eil101":   101,
    "lin105":   105,
    "pr107":    107,
    "pr124":    124,
    "bier127":  127,
    "pr136":    136,
    "pr144":    144,
    "ch150":    150,
    "kroA150":  150,
    "kroB150":  150,
    "pr152":    152,
    "u159":     159,
    "rat195":   195,
    "d198":     198,
    "kroA200":  200,
    "kroB200":  200,
    "ts225":    225,
    "tsp225":   225,
    "pr226":    226,
    "gil262":   262,
    "pr264":    264,
    "a280":     280,
    "pr299":    299,
    "lin318":   318,
    "linhp318": 318,
    "rd400":    400,
    "fl417":    417,
    "pr439":    439,
    "pcb442":   442,
    "d493":     493,
    "u574":     574,
    "rat575":   575,
    "p654":     654,
    "d657":     657,
    "u724":     724,
    "rat783":   783,
    "pr1002":   1002,
    "u1060":    1060,
    "vm1084":   1084,
    "pcb1173":  1173,
    "d1291":    1291,
    "rl1304":   1304,
    "rl1323":   1323,
    "nrw1379":  1379,
    "fl1400":   1400,
    "u1432":    1432,
    "fl1577":   1577,
    "d1655":    1655,
    "vm1748":   1748,
    "u1817":    1817,
    "rl1889":   1889,
    "d2103":    2103,
    "u2152":    2152,
    "u2319":    2319,
    "pr2392":   2392,
    "pcb3038":  3038,
    "fl3795":   3795,
    "fnl4461":  4461,
    "rl5915":   5915,
    "rl5934":   5934,
    "rl11849":  11849,
    "usa13509": 13509,
    "brd14051": 14051,
    "d15112":   15112,
    "d18512":   18512
}

solver = sys.argv[1]
process = subprocess.Popen("./main",
                           stdin=subprocess.PIPE,
                           stdout=subprocess.PIPE,
                           text=True,
                           bufsize=1,
                           universal_newlines=True)

def send(text):
    process.stdin.write(text + '\n')
    process.stdin.flush()

def finish():
    process.stdin.close()
    output_lines = []
    while True:
        line = process.stdout.readline()
        if not line:
            break
        print(line,end="")
        output_lines.append(line)
    process.stdout.close()
    process.wait()

    dt = datetime.now()
    dtstr = dt.strftime("%Y-%m-%dT%H_%M_%S")
    with open("experiments/experiments_stdout_" + dtstr + ".txt", "w") as file:
        file.write("".join(output_lines))


send("solver " + solver)
send("iterations " + str(max_iterations))
send("executions " + str(executions))
for problem in selected_problems:
    send("problem " + problem)
    cutoff = lower_bounds[problem] * cutoff_percent(dimensions[problem])
    send("cutoff " + str(cutoff))
    send("run")
finish()

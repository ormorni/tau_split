# The Tau-Splitting algorithm

The Tau-Splitting algorithm is an algorithm for the simulation of chemical reaction networks.

The code is implemented using the Rust programming language, which can be installed from the [official website](https://rust-lang.org/).
After Rust is installed, the Tau-Split algorithm can be installed by typing `cargo install tausplit` in the command line.
Afterwards, it can be used from the command line using

```bash
tausplit {time} {input_file}
```
to simulate the chemical reaction network specified in the input file for the specified time period,
and output a TSV containing the initial and final states of the system.

## Options

```bash
-s {count}, --samples {count}
```
Samples the state of the system `count` evenly spaced times. If `time` is 40 and `count` is 4, the algorithm will sample the system at 10, 20, 30, and 40 and output the states at those times in the TSV.

```bash
--algorithm {algorithm}
```

Uses the given algorithm to simulate the system. The currently supported algorithms are `tau-split`, `tau-split6`, and `gillespie`. The default algorithm is `tau-split`, which is an optimized version of the algorithm described in the manuscript.

```bash
--seed {seed}
```

Uses the given seed for random number generation during the algorithm's run.

```bash
--count-reactions
```

Adds the total number of reactions to the TSV produced.

```bash
--cpu-time
```

Adds the total runtime of the algorithm to the TSV produced.

```bash
--no-print-state
```

Makes the program not include the state of the system in the TSV produced.

## Input format

The chemical reaction networks are specified by files in the following format:

```
A = 5
B = 3
C = 0
D = 0

A + B -> C + 3D, 3.5
```

Rows of the format `A = 5` define that a chemical species `A` exists, and that its amount at the beginning of the simulation is 5.
Rows of the format `A + B -> C + 3D, 3.5` define the chemical reaction in which `A` and `B` react to form one `C` and three `D`s, with a rate constant of 3.5.
The reaction network can be split across any number of files, allowing using the same chemical reaction network with different initial states.

Lines starting with `#` are treated as comments, and are not parsed.
# Half-approximate Graph Matching in UPC++

This folder contains codes for implementing Half-approximate
matching using [UPC++](https://upcxx.lbl.gov).

Please review the following paper for the serial
algorithm, based on Manne-Bisseling:     
http://www.staff.science.uu.nl/~bisse101/Articles/CP75.pdf

This is based on our MPI-based half approximate graph matching 
implementation, discussed in paper:    
https://ieeexplore.ieee.org/abstract/document/8820975


## Compile

Just invoking `make` should build the program, without any
changes made to the Makefile, provided `upcxx` and MPI are in the `$PATH`. 
It will use whatever network backend is configured by default
for your UPC++ install, which is site-dependent.

If you have an InfiniBand network, you can override the default 
with a command like:
```
make FLAGS="-O3 -network=ibv"
```

## Running on an HPC system

On an HPC system with a proper UPC++ install, you should be able
to run with a command like:

```
mpiexec -n 2 ./match_rma -n 1000
```

Note this requires that UPC++ was built with support for MPI hybrid
applications (enabled by default on most platforms).

## Running on a laptop

On a laptop/workstation, set these prior to calling upcxx-run:

```
export UPCXX_GASNET_CONDUIT=udp
export GASNET_SPAWNFN='C'
export GASNET_CSPAWN_CMD='mpirun -np %N %C'
```

Example: 

```
upcxx-run -n 2 ./match_rma -n 1000
```

## Command line options:


Apart from using external file, it is possible to generate
in-memory a random geometric graph in a distributed fashion.
To learn more about binary file format for real-world graphs, 
see FAQs of miniVite: https://github.com/Exa-Graph/miniVite

Possible options (can be combined):
| Option          | Description                                               |
| --------------- | --------------------------------------------------------- |
| `-f <bin-file>` | Specify input binary file after this argument.            |
| `-n <vertices>` | Pass total number of vertices of the generated graph.     |
| `-l`            | Use distributed LCG for randomly choosing edges. If this option is not used, we will use C++ random number generator (using `std::default_random_engine`). |
| `-p <percent>`  | Specify percent of overall edges to be randomly generated between processes. |
| `-w`            | Use Euclidean distance as edge weight. If this option is not used, edge weights are considered as 1.0. Generate edge weight uniformly between (0,1) if Euclidean distance is not available (applicable to randomly generated edges). |
| `-r <nranks>`   | This is used to control the number of aggregators in MPI I/O and is meaningful when an input binary graph file is passed with option `-f`.  `naggr := (nranks > 1) ? (nprocs/nranks) : nranks;` 

## Contributors

* [Pagoda project](https://crd.lbl.gov/pagoda) at [LBNL](https://www.lbl.gov)
* PNNL Data Sciences



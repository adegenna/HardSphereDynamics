![unittests](https://github.com/adegenna/HardSphereDynamics/actions/workflows/cmake.yml/badge.svg)

Solver to compute 2D hard sphere collision dynamics. Parallelized with OpenMP.

Copyright 2021 by Anthony M. DeGennaro (ISC License).

**Installation**

```sh
cd [/PATH/TO/HardSphereDynamics]
mkdir build
cd build
cmake ../
make
```

**Unit Tests**

```sh
cd [/PATH/TO/HardSphereDynamics]/build
./testall
```

**Example Driver**

```sh
cd [/PATH/TO/HardSphereDynamics]/examples
./run_example.sh
./plot_example.sh
```
![](https://github.com/adegenna/HardSphereDynamics/blob/master/figs/billiards.png)
*Snapshots in time of a billiard ball rack break*

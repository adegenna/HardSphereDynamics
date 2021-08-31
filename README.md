![unittests](https://github.com/adegenna/embeddingSampler/actions/workflows/unittests.yml/badge.svg)

Library for tools related to sampling from a linear embedded subspace. Allows user to define a target distribution in an ambient space (with box constraints) of arbitrary dimension, then draw samples from an embedded linear subspace which map to that target distribution. See https://arxiv.org/abs/2001.11659 for details, as this is part of the ALEBO algorithm for high dimensional Bayesian inference.

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
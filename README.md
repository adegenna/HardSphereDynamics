![unittests](https://github.com/adegenna/HardSphereDynamics/actions/workflows/cmake.yml/badge.svg)

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

This example shows a break of a simple four-row rack of billiard balls. Here are three successive snapshots in time showing this:

<img src=https://user-images.githubusercontent.com/2964258/131427559-99fe4f70-8706-483d-9493-7a28eb01b4d3.png width="256"> <img src=https://user-images.githubusercontent.com/2964258/131427575-51f47fb3-4283-4cfb-b337-36c4de6b7738.png width="256"> <img src=https://user-images.githubusercontent.com/2964258/131427585-3b2e8cbb-96ff-47c7-abc4-d33b5cd8fe0f.png width="256">

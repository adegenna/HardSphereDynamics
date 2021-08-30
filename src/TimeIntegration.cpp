#include "TimeIntegration.h"
#include "LagrangianState.h"
#include "Inputfile.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <istream>
#include <random>
#include <cmath>
#include <boost/program_options.hpp>
#include <iostream>
#include <iterator>
#include <omp.h>
#include <sys/time.h>

using namespace std;
using namespace Eigen;
namespace po = boost::program_options;

TimeIntegration::TimeIntegration(Options& options, ParticlePhysics& physics, LagrangianState& state)
  : options_(options), physics_(&physics), state_(&state)
{
}

TimeIntegration::TimeIntegration(Options& options, ParticlePhysics& physics, LagrangianState& state, LagrangianState& boundary)
  : options_(options), physics_(&physics), state_(&state), boundary_(&boundary)
{
}

TimeIntegration::~TimeIntegration() {
  
}

void TimeIntegration::euler() {
  MatrixXd DUV, DXY;
  for (int i=0; i<options_.tsteps; i++) {
    // Calculate particle interactions: x_tt = RHS
    const MatrixXd& rhs = physics_->RHS();
    DUV = rhs*options_.dt;
    state_->incrementUV(DUV);
    const MatrixXd& UV = state_->getUV();
    DXY = UV*options_.dt;
    state_->updateXY(DXY);
    if ( (i+1)%options_.tsave == 0) {
      state_->writeXY(options_.outputfile+"_"+std::to_string(i+1)+".csv");
    }

  }

  
}

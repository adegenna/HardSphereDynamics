#include <boost/program_options.hpp>
#include "physics_test.h"
#include "../src/Inputfile.hpp"
#include "../src/LagrangianState.h"
#include "../src/ParticlePhysics.h"
#include "../src/TimeIntegration.h"
#include <sys/time.h>

using namespace Eigen;

TEST_F(PhysicsTest, testEulerIntegrator) {
  
  Options options;
  options.inputfile  = std::string(SRCDIR)+"tests/testinput.csv";
  options.outputfile = "final.csv";
  options.dt         = 0.1;
  options.tsteps     = 10;
  options.tsave      = 20;
  
  // Setup
  LagrangianState simulation(load_csv<MatrixXd>(options.inputfile));
  ParticlePhysics physics(options, simulation);
  TimeIntegration integrator(options, physics, simulation);
  
  // Solve
  integrator.euler();

  // Output
  simulation.writeXY(options.outputfile);
  
  // Read the output
  MatrixXd out;
  out   = load_csv<MatrixXd>(options.outputfile);
  
  ASSERT_TRUE(out.isApprox(Xfinal_));
  
}

TEST_F(PhysicsTest, testBilliards) {

  Options options;
  options.inputfile  = std::string(SRCDIR)+"tests/billiards.csv";
  options.outputfile = "billiardsfinal_BruteForce";
  options.dt         = 0.001;
  options.tsteps     = 10000;
  options.tsave      = 500;

  Options options_rtree;
  options_rtree.inputfile  = std::string(SRCDIR)+"tests/billiards.csv";
  options_rtree.outputfile = "billiardsfinal_rtree";
  options_rtree.dt         = 0.001;
  options_rtree.tsteps     = 10000;
  options_rtree.tsave      = 500;
  
  // Setup models
  LagrangianState simulation(load_csv<MatrixXd>(options.inputfile));
  ParticlePhysics physics(options, simulation, false);
  TimeIntegration integrator(options, physics, simulation);
  
  LagrangianState simulation_rtree(load_csv<MatrixXd>(options_rtree.inputfile));
  ParticlePhysics physics_rtree(options_rtree, simulation_rtree, true);
  TimeIntegration integrator_rtree(options_rtree, physics_rtree, simulation_rtree);
  
  // Solve
  struct timeval start, end;
  gettimeofday(&start, NULL);
  integrator.euler();
  gettimeofday(&end, NULL);
  double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + 
		  end.tv_usec - start.tv_usec) / 1.e6;
  std::cout << "Brute force calculation: " << delta << " s" << std::endl;  
  
  gettimeofday(&start, NULL);
  integrator_rtree.euler();
  gettimeofday(&end, NULL);
  delta = ((end.tv_sec  - start.tv_sec) * 1000000u + 
		  end.tv_usec - start.tv_usec) / 1.e6;  
  std::cout << "Rtree calculation: " << delta << " s" << std::endl;
  
  // Output
  simulation.writeXY(options.outputfile);
  simulation_rtree.writeXY(options_rtree.outputfile);
  
  // Read the output
  MatrixXd out, out2;
  out   = load_csv<MatrixXd>(options.outputfile);
  out2  = load_csv<MatrixXd>(options_rtree.outputfile);

  EXPECT_LT( (out-out2).norm(), 1.e-1);
  EXPECT_LT( (out-XYfinalBilliards_).norm(), 1.e-1);
  
}

TEST(ContactForceModel, OverlappingStationaryParticlesHaveZeroForce)
{
  Eigen::MatrixXd state(2,5);
  // stationary particle at 0,0 with radius 1
  state(0,0) = 0;
  state(0,1) = 0;
  state(0,2) = 0;
  state(0,3) = 0;
  state(0,4) = 1;
  // stationary particle at 1,0 with radius 1
  state(1,0) = 1;
  state(1,1) = 0;
  state(1,2) = 0;
  state(1,3) = 0;
  state(1,4) = 1;

  BruteForceContactForceModel brute_force_contact;
  RTreeContactForceModel rtree_contact;
  Eigen::MatrixXd forces_brute_force, forces_rtree;
  
  brute_force_contact.particleContact(LagrangianState(state), forces_brute_force);
  rtree_contact.particleContact(LagrangianState(state), forces_rtree);

  EXPECT_TRUE(forces_brute_force.isApprox(forces_rtree));
  
  Eigen::MatrixXd zero2x2 = Eigen::MatrixXd::Zero(2,2);
  EXPECT_TRUE(forces_brute_force.isApprox(zero2x2)) << "brute force should be zero.";
  EXPECT_TRUE(forces_rtree.isApprox(zero2x2)) << "rtree force should be zero.";
}


TEST(ContactForceModel, OverlappingInX)
{
  Eigen::MatrixXd state(2,5);
  // particle at 0,0 with radius 0.51 moving at (1,0)
  state(0,0) = 0;
  state(0,1) = 0;
  state(0,2) = 1;
  state(0,3) = 0;
  state(0,4) = 0.51;
  // particle at 1,0 with radius 0.51 moving at (-1,0)
  state(1,0) = 1;
  state(1,1) = 0;
  state(1,2) = -1;
  state(1,3) = 0;
  state(1,4) = 0.51;

  BruteForceContactForceModel brute_force_contact;
  RTreeContactForceModel rtree_contact;
  Eigen::MatrixXd forces_brute_force, forces_rtree;

  brute_force_contact.particleContact(LagrangianState(state), forces_brute_force);
  rtree_contact.particleContact(LagrangianState(state), forces_rtree);

  EXPECT_TRUE(forces_brute_force.isApprox(forces_rtree));

  Eigen::MatrixXd zero2x2 = Eigen::MatrixXd::Zero(2,2);
  EXPECT_FALSE(forces_brute_force.isApprox(zero2x2)) << "brute force should not be zero.";
  EXPECT_FALSE(forces_rtree.isApprox(zero2x2)) << "rtree force should not be zero.";

  Eigen::MatrixXd expected(2,2);
  expected(0,0) = -286.622*3.77;
  expected(1,0) = 286.622*3.77;
  expected(0,1) = 0;
  expected(1,1) = 0;

  EXPECT_TRUE(forces_brute_force.isApprox(expected, 0.1));
  EXPECT_TRUE(forces_rtree.isApprox(expected, 0.1));
}
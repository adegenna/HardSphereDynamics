#ifndef PHYSICS_TEST_H_
#define PHYSICS_TEST_H_

#include "gtest/gtest.h"
#include "math.h"
#include <Eigen/Dense>

class PhysicsTest: public ::testing::Test {
 protected:
  virtual void SetUp() {
    Xfinal_.resize(4,2);
    Xfinal_ << 11.0, 11.0,
      -11.0, 11.0,
      -11.0, -11.0,
      11.0, -11.0;
    XYfinalBilliards_.resize(11,2);
    XYfinalBilliards_ << 0, 0.755669 ,
      -1.9791, 3.2347 , 
      1.9791, 3.2347 , 
      -3.39979, 6.07807 , 
      0, 4.73963 , 
      3.39979, 6.07807 , 
      -4.37855, 9.3162 , 
      -1.24109, 8.88838 , 
      1.24109, 8.88838 , 
      4.37855, 9.3162 , 
      0, -2.56848;
  }

  Eigen::MatrixXd Xfinal_;
  Eigen::MatrixXd XYfinalBilliards_;
  
};

#endif

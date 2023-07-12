#include "common.hpp"

TEST(UnitTest, ortho2D) {
  ortho<2,float> R1 = rotation_matrix(0.25f);
  EXPECT_NEAR(det(as_mat(R1)), 1.0f, epsf);

  ortho<2,double> R2 = rotation_matrix(0.25);
  EXPECT_NEAR(det(as_mat(R2)), 1.0, eps);
}

TEST(UnitTest, ortho3D) {
  ortho<3,float> R1 = rotation_matrix(vec3f{0.25f, 0.1f, 0.4f});
  EXPECT_NEAR(det(as_mat(R1)), 1.0f, 2.0 * epsf);

  ortho<3,double> R2 = rotation_matrix(vec3{0.25, 0.1, 0.4});
  EXPECT_NEAR(det(as_mat(R2)), 1.0, eps);
  EXPECT_NEAR(det(as_mat(transpose(R2))), 1.0, eps);
}
#include "common.hpp"

#include "fm/operations/dot.hpp"
#include "fm/operations/linear_solve.hpp"

double tolerance = 1.0e-14;

diag<2, double> D2{{4.0, 2.0}};
diag<3, double> D3{{2.0, 4.0, 8.0}};

iso<2, double> I2{2.0};
iso<3, double> I3{4.0};

sym<2, double> S2{{
  5.0, 3.0,
       2.0
}};

sym<3, double> S3{{{
  -2.0,  1.0,  0.0,
        -2.0,  1.0,
              -2.0
}}};

mat<2,2,double> M2{{
  {5, 3},
  {3, 2}
}};

mat<2,2,double> invM2{{
  {2, -3}, 
  {-3, 5}
}};

mat<3,3,double> M3{{
  {-2.0,  1.0,  0.0},
  { 1.0, -2.0,  1.0},
  { 0.0,  1.0, -2.0}
}};

mat<3,3,double> invM3{{
  {-0.75, -0.5, -0.25}, 
  {-0.50, -1.0, -0.50}, 
  {-0.25, -0.5, -0.75}
}};

vec2 b2 = {1.0, 2.0};
vec3 b3 = {1.0, 2.0, 3.0};

TEST(UnitTest, linear_solve_mat2) {
  vec2 x2 = linear_solve(M2, b2);
  vec2 expected = dot(invM2, b2);
  for (int i = 0; i < 2; i++) {
    EXPECT_NEAR(x2[i], expected[i], tolerance);
  }
}

TEST(UnitTest, linear_solve_mat3) {
  vec3 x3 = linear_solve(M3, b3);
  vec3 expected = dot(invM3, b3);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(x3[i], expected[i], tolerance);
  }
}

TEST(UnitTest, linear_solve_sym2) {
  vec2 x2 = linear_solve(S2, b2);
  vec2 expected = dot(invM2, b2);
  for (int i = 0; i < 2; i++) {
    EXPECT_NEAR(x2[i], expected[i], tolerance);
  }
}

TEST(UnitTest, linear_solve_sym3) {
  vec3 x3 = linear_solve(S3, b3);
  vec3 expected = dot(invM3, b3);
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(x3[i], expected[i], tolerance);
  }
}

TEST(UnitTest, linear_solve_diag2) {
  vec2 x2 = linear_solve(D2, b2);
  vec2 expected = {b2[0] / D2.data[0], b2[1] / D2.data[1]};
  for (int i = 0; i < 2; i++) {
    EXPECT_NEAR(x2[i], expected[i], tolerance);
  }
}

TEST(UnitTest, linear_solve_diag3) {
  vec3 x3 = linear_solve(D3, b3);
  vec3 expected = {b3[0] / D3.data[0], b3[1] / D3.data[1], b3[2] / D3.data[2]};
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(x3[i], expected[i], tolerance);
  }
}

TEST(UnitTest, linear_solve_iso2) {
  vec2 x2 = linear_solve(I2, b2);
  vec2 expected = {b2[0] / I2.data, b2[1] / I2.data};
  for (int i = 0; i < 2; i++) {
    EXPECT_NEAR(x2[i], expected[i], tolerance);
  }
}

TEST(UnitTest, linear_solve_iso3) {
  vec3 x3 = linear_solve(I3, b3);
  vec3 expected = {b3[0] / I3.data, b3[1] / I3.data, b3[2] / I3.data};
  for (int i = 0; i < 3; i++) {
    EXPECT_NEAR(x3[i], expected[i], tolerance);
  }
}

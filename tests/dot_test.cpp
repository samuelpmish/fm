#include "common.hpp"

using namespace fm;

iso2 I2{4.0};
diag2 D2{{2.0, 3.0}};
skew2 A2{1.5};
rot2 R2 = RotationMatrix(1.0);
sym2 S2{{1.0, 2.0, 3.0}};
mat2 M2{{{1.0, 2.0}, {3.0, 4.0}}};

iso3 I3{2.0};
diag3 D3{{2.0, 3.0, 4.0}};
skew3 A3{{-1.5, 1.0, -2.0}};
rot3 R3 = RotationMatrix(vec3{1.0, 0.1, 0.5});
sym3 S3{{{1.0, 2.0, 3.0,
               4.0, 5.0,
                    6.0}}};
mat3 M3{{{1.0, 2.0, 3.0}, 
         {4.0, 5.0, 6.0}, 
         {7.0, 8.0, 10.0}}};

TEST(dot, iso2_iso2) {
  EXPECT_NEAR(norm(dot(I2, I2) - iso2{16.0}), 0.0, 1.0e-15);
}

TEST(dot, iso2_diag2) {
  EXPECT_NEAR(norm(dot(I2, D2) - diag2{{8.0, 12.0}}), 0.0, 1.0e-15);
}

TEST(dot, iso2_skew2) {
  EXPECT_NEAR(norm(dot(I2, A2) - skew2{6}), 0.0, 1.0e-15);
}

TEST(dot, iso2_rot2) {
  EXPECT_NEAR(norm(dot(I2, R2) - 4.0 * mat2(R2)), 0.0, 1.0e-15);
}

TEST(dot, iso2_sym2) {
  EXPECT_NEAR(norm(dot(I2, S2) - 4.0 * mat2(S2)), 0.0, 1.0e-15);
}

TEST(dot, iso2_mat2) {
  EXPECT_NEAR(norm(dot(I2, M2) - 4.0 * M2), 0.0, 1.0e-15);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

TEST(dot, diag2_iso2) {
  EXPECT_NEAR(norm(dot(D2, I2) - diag2{{8.0, 12.0}}), 0.0, 1.0e-15);
}

TEST(dot, diag2_diag2) {
  EXPECT_NEAR(norm(dot(D2, D2) - diag2{{4.0, 9.0}}), 0.0, 1.0e-15);
}

TEST(dot, diag2_skew2) {
  EXPECT_NEAR(norm(dot(D2, A2) - mat2{{{0.0, -3.0}, {4.5, 0.0}}}), 0.0, 1.0e-15);
}

TEST(dot, diag2_rot2) {
  double c1 = cos(1.0);
  double s1 = sin(1.0);
  EXPECT_NEAR(norm(dot(D2, R2) - mat2{{{2*c1, -2*s1},{3*s1, 3*c1}}}), 0.0, 1.0e-15);
}

TEST(dot, diag2_sym2) {
  EXPECT_NEAR(norm(dot(D2, S2) - mat2{{{2.0, 4.0},{6.0, 9.0}}}), 0.0, 1.0e-15);
}

TEST(dot, diag2_mat2) {
  EXPECT_NEAR(norm(dot(D2, M2) - mat2{{{2.0, 4.0}, {9.0, 12.0}}}), 0.0, 1.0e-15);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

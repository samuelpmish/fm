#include "common.hpp"

#include "types/matrix.hpp"
#include "operations/print.hpp"
#include "operations/inverse.hpp"

using namespace fm;

iso2 I2{4.0};
diag2 D2{2.0, 3.0};
skew2 A2{-1.5};
rot2 R2 = RotationMatrix(1.0);
sym2 S2{1.0, 2.0, 3.0};
mat2 M2{{{1.0, 2.0}, {3.0, 4.0}}};

TEST(inverse, iso2) {
  mat2 expected = {{{0.25, 0.0}, {0.0, 0.25}}};
  EXPECT_NEAR(norm(inv(I2) - expected), 0.0, 1.0e-15);
  EXPECT_NEAR(norm(inv(I2) - iso2{0.25}), 0.0, 1.0e-15);
}

TEST(inverse, diag2) {
  mat2 expected = {{{0.5, 0.0}, {0.0, 1.0 / 3.0}}};
  EXPECT_NEAR(norm(inv(D2) - expected), 0.0, 1.0e-15);
  EXPECT_NEAR(norm(inv(D2) - diag2{0.5, 1.0 / 3.0}), 0.0, 1.0e-15);
}

TEST(inverse, skew2) {
  mat2 expected = {{{0.0, 2.0 / 3.0}, {-2.0 / 3.0, 0.0}}};
  EXPECT_NEAR(norm(inv(A2) - expected), 0.0, 1.0e-15);
  EXPECT_NEAR(norm(inv(A2) - skew2{2.0 / 3.0}), 0.0, 1.0e-15);
}

TEST(inverse, rot2) {
  mat2 expected = RotationMatrix(-1.0);
  EXPECT_NEAR(norm(inv(R2) - expected), 0.0, 1.0e-15);
}

TEST(inverse, sym2) {
  mat2 expected{{{-3, 2}, {2, -1}}};
  EXPECT_NEAR(norm(inv(S2) - expected), 0.0, 1.0e-15);
}

TEST(inverse, mat2) {
  mat2 expected = {{{-2, 1}, {1.5, -0.5}}};
  EXPECT_NEAR(norm(inv(M2) - expected), 0.0, 1.0e-15);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

iso3 I3{2.0};
diag3 D3{2.0, 3.0, 4.0};
skew3 A3{-1.5, 1.0, -2.0};
rot3 R3 = RotationMatrix(vec3{1.0, 0.1, 0.5});
sym3 S3{1.0, 2.0, 3.0,
             4.0, 5.0,
                  6.0};
mat3 M3{{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 10.0}}};

TEST(inverse, iso3) {
  mat3 expected = {{{0.5, 0.0, 0.0}, {0.0, 0.5, 0.0}, {0.0, 0.0, 0.5}}};
  EXPECT_NEAR(norm(inv(I3) - expected), 0.0, 1.0e-15);
  EXPECT_NEAR(norm(inv(I3) - iso3{0.5}), 0.0, 1.0e-15);
}

TEST(inverse, diag3) {
  mat3 expected = {{{0.5, 0.0, 0.0}, {0.0, 1.0 / 3.0, 0.0}, {0.0, 0.0, 0.25}}};
  EXPECT_NEAR(norm(inv(D3) - expected), 0.0, 1.0e-15);
  EXPECT_NEAR(norm(inv(D3) - diag3{0.5, 1.0 / 3.0, 0.25}), 0.0, 1.0e-15);
}

// this test should fail at compile time, 
// since all skew3 matrices are singular
//TEST(inverse, skew3) { inv(A3); }

TEST(inverse, rot3) {
  mat3 expected = RotationMatrix(vec3{-1.0, -0.1, -0.5});
  EXPECT_NEAR(norm(inv(R3) - expected), 0.0, 1.0e-15);
}

TEST(inverse, sym3) {
  mat3 expected{{{1, -3, 2}, {-3, 3, -1}, {2, -1, 0}}};
  EXPECT_NEAR(norm(inv(S3) - expected), 0.0, 1.0e-15);
}

TEST(inverse, mat3) {
  mat3 expected = {{{-(2.0/3.0), -(4.0/3.0), 1.0}, 
                    {-(2.0/3.0), 11.0/3.0, -2.0}, 
                    {1.0, -2.0, 1.0}}};
  EXPECT_NEAR(norm(inv(M3) - expected), 0.0, 1.0e-15);
}
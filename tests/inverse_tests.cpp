#include "common.hpp"

#include "types/matrix.hpp"
#include "operations/print.hpp"
#include "operations/inverse.hpp"

namespace fm {

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

}
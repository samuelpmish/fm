#include "common.hpp"

TEST(invariants, tr_mat2) {
  EXPECT_NEAR(tr(mat2{{{1.0, 2.0}, {3.0, 4.0}}}), 5.0, 1.0e-15);
}
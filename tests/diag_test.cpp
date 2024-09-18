#include "common.hpp"

using namespace fm;

TEST(UnitTest, diag2D) {
  diag<2> I1{1.0, 2.0};
  diag<2> I2{0.25, 0.5};

  EXPECT_EQ((I1 + I2).data, vec({1.25, 2.5}));
  EXPECT_EQ((I1 - I2).data, vec({0.75, 1.5}));
  EXPECT_EQ((I1 * I2).data, vec({0.25, 1.0}));
  EXPECT_EQ((I1 / I2).data, vec({4.0, 4.0}));
  EXPECT_EQ(dot(I1, I2).data, vec({0.25, 1.0}));
  EXPECT_EQ((I1 * inv(I2)).data, vec({4.0, 4.0}));
  EXPECT_EQ(det(I1 * inv(I2)), 16.0);
}

TEST(UnitTest, diag3D) {
  diag<3> I1{1.0, 2.0, 4.0};
  diag<3> I2{0.25, 0.5, 1.0};

  EXPECT_EQ((I1 + I2).data, vec({1.25, 2.5, 5.0}));
  EXPECT_EQ((I1 - I2).data, vec({0.75, 1.5, 3.0}));
  EXPECT_EQ((I1 * I2).data, vec({0.25, 1.0, 4.0}));
  EXPECT_EQ((I1 / I2).data, vec({4.0, 4.0, 4.0}));
  EXPECT_EQ(dot(I1, I2).data, vec({0.25, 1.0, 4.0}));
  EXPECT_EQ((I1 * inv(I2)).data, vec({4.0, 4.0, 4.0}));
  EXPECT_EQ(det(I1 * inv(I2)), 64.0);
}

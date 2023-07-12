#include "common.hpp"

TEST(UnitTest, diag2D) {
  diag<2, float> I1{1.0f, 2.0f};
  diag<2, float> I2{0.25f, 0.5f};

  EXPECT_EQ((I1 + I2).data, vec({1.25f, 2.5f}));
  EXPECT_EQ((I1 - I2).data, vec({0.75f, 1.5f}));
  EXPECT_EQ((I1 * I2).data, vec({0.25f, 1.0f}));
  EXPECT_EQ((I1 / I2).data, vec({4.0f, 4.0f}));
  EXPECT_EQ(dot(I1, I2).data, vec({0.25f, 1.0f}));
  EXPECT_EQ((I1 * inv(I2)).data, vec({4.0f, 4.0f}));
  EXPECT_EQ(det(I1 * inv(I2)), 16.0f);
}

TEST(UnitTest, diag3D) {
  diag<3, float> I1{1.0f, 2.0f, 4.0f};
  diag<3, float> I2{0.25f, 0.5f, 1.0f};

  EXPECT_EQ((I1 + I2).data, vec({1.25f, 2.5f, 5.0f}));
  EXPECT_EQ((I1 - I2).data, vec({0.75f, 1.5f, 3.0f}));
  EXPECT_EQ((I1 * I2).data, vec({0.25f, 1.0f, 4.0f}));
  EXPECT_EQ((I1 / I2).data, vec({4.0f, 4.0f, 4.0f}));
  EXPECT_EQ(dot(I1, I2).data, vec({0.25f, 1.0f, 4.0f}));
  EXPECT_EQ((I1 * inv(I2)).data, vec({4.0f, 4.0f, 4.0f}));
  EXPECT_EQ(det(I1 * inv(I2)), 64.0f);
}

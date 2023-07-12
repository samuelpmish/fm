#include "common.hpp"

TEST(UnitTest, ISO2D) {
  iso<2, float> I1{2.0f};
  iso<2, float> I2{0.25f};

  EXPECT_EQ((I1 + I2).data, 2.25f);
  EXPECT_EQ((I1 - I2).data, 1.75f);
  EXPECT_EQ((I1 * I2).data, 0.5f);
  EXPECT_EQ((I1 / I2).data, 8.0f);
  EXPECT_EQ(dot(I1, I2).data, 0.5f);
  EXPECT_EQ((I1 * inv(I2)).data, 8.0f);
  EXPECT_EQ(det(I1 * inv(I2)), 64.0f);
}

TEST(UnitTest, ISO3D) {
  iso<3, float> I1{2.0f};
  iso<3, float> I2{0.25f};

  EXPECT_EQ((I1 + I2).data, 2.25f);
  EXPECT_EQ((I1 - I2).data, 1.75f);
  EXPECT_EQ((I1 * I2).data, 0.5f);
  EXPECT_EQ((I1 / I2).data, 8.0f);
  EXPECT_EQ(dot(I1, I2).data, 0.5f);
  EXPECT_EQ((I1 * inv(I2)).data, 8.0f);
  EXPECT_EQ(det(I1 * inv(I2)), 512.0f);
}

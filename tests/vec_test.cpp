#include "common.hpp"

TEST(UnitTest, vec2D) {
  vec2f u = normalize(vec2f{0.25f, 0.1f});
  EXPECT_NEAR(u[0], 0.9284766908852594, epsf);
  EXPECT_NEAR(u[1], 0.3713906763541038, epsf);
}

TEST(UnitTest, vec3D) {
  vec3f u = normalize(vec3f{0.25f, 0.1f, 0.4f});
  EXPECT_NEAR(u[0], 0.5184758473652126, epsf);
  EXPECT_NEAR(u[1], 0.2073903389460851, epsf);
  EXPECT_NEAR(u[2], 0.8295613557843402, epsf);
}

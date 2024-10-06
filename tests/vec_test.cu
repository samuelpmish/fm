#include "common.hpp"

__global__ void test_kernel() {

  pass = true;

  {
    vec2f u = normalize(vec2f{0.25f, 0.1f});
    CUDA_EXPECT_NEAR(u[0], 0.9284766908852594, epsf);
    CUDA_EXPECT_NEAR(u[1], 0.3713906763541038, epsf);
  }

  {
    vec3f u = normalize(vec3f{0.25f, 0.1f, 0.4f});
    CUDA_EXPECT_NEAR(u[0], 0.5184758473652126, epsf);
    CUDA_EXPECT_NEAR(u[1], 0.2073903389460851, epsf);
    CUDA_EXPECT_NEAR(u[2], 0.8295613557843402, epsf);
  }

  if (!pass) { asm("trap;"); }

} 

TEST(cuda, vec_tests){
  test_kernel<<<1,1>>>();
  EXPECT_FALSE(cudaDeviceSynchronize());
}
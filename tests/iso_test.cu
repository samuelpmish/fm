#include "common.hpp"

__global__ void test_kernel() {

  pass = true;

  {
    iso<2, float> I1{2.0f};
    iso<2, float> I2{0.25f};

    CUDA_EXPECT_EQ(dot(I1, I2).data, 0.5f);
    CUDA_EXPECT_EQ(dot(I1, inv(I2)).data, 8.0f);
    CUDA_EXPECT_EQ(det(dot(I1, inv(I2))), 64.0f);
  }

  {
    iso<3, float> I1{2.0f};
    iso<3, float> I2{0.25f};

    CUDA_EXPECT_EQ(dot(I1, I2).data, 0.5f);
    CUDA_EXPECT_EQ(dot(I1, inv(I2)).data, 8.0f);
    CUDA_EXPECT_EQ(det(dot(I1, inv(I2))), 512.0f);
  }

  if (!pass) { asm("trap;"); }

} 

TEST(cuda, iso_tests){
  test_kernel<<<1,1>>>();
  EXPECT_FALSE(cudaDeviceSynchronize());
}

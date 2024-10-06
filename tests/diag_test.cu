#include "common.hpp"

using namespace fm;

__global__ void test_kernel() {

  pass = true;

  {
    diag<2> D1{{1.0, 2.0}};
    diag<2> D2{{0.25, 0.5}};

    CUDA_EXPECT_EQ((D1 + D2).data, vec({1.25, 2.5}));
    CUDA_EXPECT_EQ((D1 - D2).data, vec({0.75, 1.5}));
    CUDA_EXPECT_EQ(dot(D1, D2).data, vec({0.25, 1.0}));
    CUDA_EXPECT_EQ(dot(D1, inv(D2)).data, vec({4.0, 4.0}));
    CUDA_EXPECT_EQ(det(dot(D1, inv(D2))), 16.0);
  }

  {
    diag<3> D1{{1.0, 2.0, 4.0}};
    diag<3> D2{{0.25, 0.5, 1.0}};

    CUDA_EXPECT_EQ((D1 + D2).data, vec({1.25, 2.5, 5.0}));
    CUDA_EXPECT_EQ((D1 - D2).data, vec({0.75, 1.5, 3.0}));
    CUDA_EXPECT_EQ(dot(D1, D2).data, vec({0.25, 1.0, 4.0}));
    CUDA_EXPECT_EQ(dot(D1, inv(D2)).data, vec({4.0, 4.0, 4.0}));
    CUDA_EXPECT_EQ(det(dot(D1, inv(D2))), 64.0);
  }

  if (!pass) { asm("trap;"); }

} 

TEST(cuda, diag_tests){
  test_kernel<<<1,1>>>();
  EXPECT_FALSE(cudaDeviceSynchronize());
}


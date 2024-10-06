#include "common.hpp"

__global__ void test_kernel() {

  pass = true;

  {
    sym<2, double> A2{{
      1.0, 2.0,
           3.0
    }};
    CUDA_EXPECT_EQ(A2.index(0, 0), 0);
    CUDA_EXPECT_EQ(A2.index(0, 1), 1);
    CUDA_EXPECT_EQ(A2.index(1, 1), 2);

    CUDA_EXPECT_EQ(A2.data[0], 1);
    CUDA_EXPECT_EQ(A2.data[1], 2);
    CUDA_EXPECT_EQ(A2.data[2], 3);
  }

  {
    sym<3, double> A3{{{
      1.0, 2.0, 3.0,
           4.0, 5.0,
                6.0
    }}};

    CUDA_EXPECT_EQ(A3.index(0, 0), 0);
    CUDA_EXPECT_EQ(A3.index(0, 1), 1);
    CUDA_EXPECT_EQ(A3.index(0, 2), 2);
    CUDA_EXPECT_EQ(A3.index(1, 1), 3);
    CUDA_EXPECT_EQ(A3.index(1, 2), 4);
    CUDA_EXPECT_EQ(A3.index(2, 1), 4);
    CUDA_EXPECT_EQ(A3.index(2, 2), 5);

    for (uint32_t i = 0; i < 6; i++) {
      CUDA_EXPECT_EQ(A3.data[i], i+1);
    }

    CUDA_EXPECT_EQ(A3(0,0) + A3(1,0) + A3(2, 0), 6);

  }

  {
    sym<6, double> A6 = {{{
      1.0, 2.0,  3.0,  4.0,  5.0,  6.0,
           7.0,  8.0,  9.0, 10.0, 11.0,
                12.0, 13.0, 14.0, 15.0,
                      16.0, 17.0, 18.0,
                            19.0, 20.0,
                                  21.0
    }}};
    CUDA_EXPECT_EQ(A6.index(0, 0),  0);
    CUDA_EXPECT_EQ(A6.index(0, 1),  1);
    CUDA_EXPECT_EQ(A6.index(0, 2),  2);
    CUDA_EXPECT_EQ(A6.index(0, 3),  3);
    CUDA_EXPECT_EQ(A6.index(0, 4),  4);
    CUDA_EXPECT_EQ(A6.index(0, 5),  5);

    CUDA_EXPECT_EQ(A6.index(1, 1),  6);
    CUDA_EXPECT_EQ(A6.index(1, 2),  7);
    CUDA_EXPECT_EQ(A6.index(1, 3),  8);
    CUDA_EXPECT_EQ(A6.index(1, 4),  9);
    CUDA_EXPECT_EQ(A6.index(1, 5), 10);

    CUDA_EXPECT_EQ(A6.index(2, 2), 11);
    CUDA_EXPECT_EQ(A6.index(2, 3), 12);
    CUDA_EXPECT_EQ(A6.index(2, 4), 13);
    CUDA_EXPECT_EQ(A6.index(2, 5), 14);

    CUDA_EXPECT_EQ(A6.index(3, 3), 15);
    CUDA_EXPECT_EQ(A6.index(3, 4), 16);
    CUDA_EXPECT_EQ(A6.index(3, 5), 17);

    CUDA_EXPECT_EQ(A6.index(4, 4), 18);
    CUDA_EXPECT_EQ(A6.index(4, 5), 19);

    CUDA_EXPECT_EQ(A6.index(5, 5), 20);

    for (uint32_t i = 0; i < 21; i++) {
      CUDA_EXPECT_EQ(A6.data[i], i+1);
    }

    CUDA_EXPECT_EQ(A6(0,0) + A6(1,0) + A6(2, 0), 6);

  }

  if (!pass) { asm("trap;"); }

} 

TEST(cuda, sym_tests){
  test_kernel<<<1,1>>>();
  EXPECT_FALSE(cudaDeviceSynchronize());
}
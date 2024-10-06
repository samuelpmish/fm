#include "common.hpp"

#include "fm/types/matrix.hpp"
#include "fm/operations/print.hpp"
#include "fm/operations/invariants.hpp"

using namespace fm;

__global__ void test_kernel() {

  pass = true;

  iso2 I2{4.0};
  diag2 D2{{2.0, 3.0}};
  skew2 A2{-1.5};
  rot2 R2 = RotationMatrix(1.0);
  sym2 S2{{1.0, 2.0, 3.0}};
  mat2 M2{{{1.0, 2.0}, {3.0, 4.0}}};

  {
    CUDA_EXPECT_NEAR(tr(I2), 8.0, 1.0e-15);
    CUDA_EXPECT_NEAR(det(I2), 16.0, 1.0e-15);
  }

  {
    CUDA_EXPECT_NEAR(tr(D2), 5.0, 1.0e-15);
    CUDA_EXPECT_NEAR(det(D2), 6.0, 1.0e-15);
  }

  {
    CUDA_EXPECT_NEAR(tr(A2), 0.0, 1.0e-15);
    CUDA_EXPECT_NEAR(det(A2), 2.25, 1.0e-15);
  }

  {
    CUDA_EXPECT_NEAR(tr(R2), tr(mat2(R2)), 1.0e-15);
    CUDA_EXPECT_NEAR(tr(R2), 2.0 * cos(1.0), 1.0e-15);
    CUDA_EXPECT_NEAR(det(R2), 1.0, 1.0e-15);
    CUDA_EXPECT_NEAR(det(R2), det(mat2(R2)), 1.0e-15);
  }

  {
    CUDA_EXPECT_NEAR(tr(S2), 4.0, 1.0e-15);
    CUDA_EXPECT_NEAR(det(S2), -1.0, 1.0e-15);
  }

  {
    CUDA_EXPECT_NEAR(tr(M2), 5.0, 1.0e-15);
    CUDA_EXPECT_NEAR(det(M2), -2.0, 1.0e-15);
  }

  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  iso3 I3{2.0};
  diag3 D3{{2.0, 3.0, 4.0}};
  skew3 A3{{-1.5, 1.0, -2.0}};
  rot3 R3 = RotationMatrix(vec3{1.0, 0.1, 0.5});
  sym3 S3{{{1.0, 2.0, 3.0,
                 4.0, 5.0,
                      6.0}}};
  mat3 M3{{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 10.0}}};

  {
    CUDA_EXPECT_NEAR(tr(I3), 6.0, 1.0e-15);
    CUDA_EXPECT_NEAR(invariant<2>(I3), -12, 1.0e-15);
    CUDA_EXPECT_NEAR(det(I3), 8.0, 1.0e-15);
  }

  {
    CUDA_EXPECT_NEAR(tr(D3), 9.0, 1.0e-15);
    CUDA_EXPECT_NEAR(invariant<1>(D3), tr(D3), 1.0e-15);
    CUDA_EXPECT_NEAR(invariant<2>(D3), -26, 1.0e-15);
    CUDA_EXPECT_NEAR(invariant<3>(D3), det(D3), 1.0e-15);
    CUDA_EXPECT_NEAR(det(D3), 24.0, 1.0e-15);
  }

  { 
    CUDA_EXPECT_NEAR(tr(A3), 0.0, 1.0e-15);
    CUDA_EXPECT_NEAR(invariant<1>(A3), tr(A3), 1.0e-15);
    CUDA_EXPECT_NEAR(invariant<2>(A3), -7.25, 1.0e-15);
    CUDA_EXPECT_NEAR(invariant<3>(A3), det(A3), 1.0e-15);
    CUDA_EXPECT_NEAR(det(A3), 0.0, 1.0e-15);
  }

  {
    CUDA_EXPECT_NEAR(tr(R3), 1.866866689763622, 1.0e-15);
    CUDA_EXPECT_NEAR(invariant<1>(R3), tr(R3), 1.0e-15);
    CUDA_EXPECT_NEAR(invariant<2>(R3), -1.866866689763622, 1.0e-15);
    CUDA_EXPECT_NEAR(invariant<3>(R3), det(R3), 1.0e-15);
    CUDA_EXPECT_NEAR(det(R3), 1.0, 1.0e-15);
  }

  {
    CUDA_EXPECT_NEAR(tr(S3), 11, 1.0e-15);
    CUDA_EXPECT_NEAR(invariant<1>(S3), tr(S3), 1.0e-15);
    CUDA_EXPECT_NEAR(invariant<2>(S3), 4, 1.0e-15);
    CUDA_EXPECT_NEAR(invariant<3>(S3), det(S3), 1.0e-15);
    CUDA_EXPECT_NEAR(det(S3), -1.0, 1.0e-15);
  }

  {
    CUDA_EXPECT_NEAR(tr(M3), 16, 1.0e-15);
    CUDA_EXPECT_NEAR(invariant<1>(M3), tr(M3), 1.0e-15);
    CUDA_EXPECT_NEAR(invariant<2>(M3), 12, 1.0e-15);
    CUDA_EXPECT_NEAR(invariant<3>(M3), det(M3), 1.0e-15);
    CUDA_EXPECT_NEAR(det(M3), -3.0, 1.0e-15);
  }

  if (!pass) { asm("trap;"); }

} 

TEST(cuda, invariant_tests){
  test_kernel<<<1,1>>>();
  EXPECT_FALSE(cudaDeviceSynchronize());
}
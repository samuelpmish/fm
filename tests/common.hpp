#pragma once

#include <gtest/gtest.h>

#include "fm.hpp"
#include "operations/print.hpp"

static constexpr float eps = 1.0e-15;
static constexpr float epsf = 1.0e-7f;

template < uint32_t m, uint32_t n, typename T >
void compare(const mat<m,n,T> & A, const mat<m,n,T> & B, double tolerance = eps) {
  for (uint32_t i = 0; i < m; i++) {
    for (uint32_t j = 0; j < n; j++) {
      EXPECT_NEAR(A[i][j], B[i][j], tolerance * std::max(0.1, abs(A[i][j]) + abs(B[i][j])));
    }
  }
}

template < uint32_t n, typename callable >
auto make_vec(callable f) {
  using T = decltype(f(int{}));
  vec<n,T> output;
  for (uint32_t i = 0; i < n; i++) {
    output[i] = f(i);
  }
  return output;
}

template < uint32_t m, uint32_t n, typename callable >
auto make_mat(callable f) {
  using T = decltype(f(int{}, int{}));
  mat<m, n,T> output;
  for (uint32_t i = 0; i < m; i++) {
    for (uint32_t j = 0; j < n; j++) {
      output[i][j] = f(i,j);
    }
  }
  return output;
}
#pragma once

#include <random>

template < typename T = double >
T random_real(T a = -1.0, T b = 1.0) {
  static std::default_random_engine generator;
  static std::uniform_real_distribution<T> distribution(0.0, 1.0);
  return a + (b - a) * distribution(generator);
}

template < uint32_t n, typename T = double >
vec<n,T> random_vec() {
  vec<n,T> output;
  for (uint32_t i = 0; i < n; i++) {
    output[i] = random_real<T>();
  }
  return output;
}

template < uint32_t m, uint32_t n, typename T = double >
mat<m,n,T> random_mat() {
  mat<m,n,T> output;
  for (uint32_t i = 0; i < n; i++) {
    output[i] = random_vec<m,T>();
  }
  return output;
}
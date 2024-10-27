#pragma once

#include <iostream>

#include "fm/macros.hpp"
#include "fm/types/matrix.hpp"

namespace fm {

template < typename T, uint32_t n >
std::ostream & operator<<(std::ostream & out, vec<n,T> v) {
  out << '{' << v(0);
  for (int i = 1; i < n; i++) {
    out << ", " << v(i);
  }
  out << '}';
  return out;
}

template < uint32_t m, uint32_t n, typename T >
std::ostream& operator<<(std::ostream& out, const mat<m,n,T> & A) {
  out << '{';
  for (uint32_t i = 0; i < m; i++) {
    out << A[i];
    if (i != n-1) out << ',';
  }
  out << '}';

  return out;
}

////////////////////////////////////////////////////////////////////////////////

__host__ __device__ inline void print(uint32_t u) { printf("%d", u); }
__host__ __device__ inline void print(float f) { printf("%f", f); }
__host__ __device__ inline void print(uint64_t u) { printf("%lu", u); }
__host__ __device__ inline void print(double f) { printf("%f", f); }

template < typename T, uint32_t n >
__host__ __device__ void print(const vec<n,T> & v) {
  printf("{");
  for (int i = 0; i < n; i++) {
    print(v[i]);
    if (i != n-1) printf(", ");
  }
  printf("}");
}

template < uint32_t m, uint32_t n, typename T >
__host__ __device__ void print(const mat<m,n,T> & A) {
  printf("{");
  for (uint32_t i = 0; i < m; i++) {
    print(A[i]);
    if (i != m-1) printf(", ");
  }
  printf("}");
}

}
#pragma once

#include "types/vec.hpp"
#include "types/mat.hpp"

template < uint32_t m, uint32_t n, typename S, typename T >
auto outer(const vec<m,S> & u, const vec<n,T> & v) {
  mat<m,n,decltype(S{} * T{})> uvT;
  for (uint32_t i = 0; i < m; i++) uvT[i] = u[i] * v;
  return uvT;
}
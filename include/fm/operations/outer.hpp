#pragma once

#include "fm/types/vec.hpp"
#include "fm/types/matrix.hpp"

namespace fm {

template < uint32_t m, uint32_t n, typename S, typename T >
auto outer(const vec<m,S> & u, const vec<n,T> & v) {
  mat<m,n,decltype(S{} * T{})> uvT;
  for (uint32_t i = 0; i < m; i++) uvT[i] = u[i] * v;
  return uvT;
}

}
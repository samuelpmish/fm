#pragma once

#include <iostream>

#include "fm/types/AABB.hpp"
#include "fm/types/vec.hpp"
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

template < typename T, int n >
std::ostream & operator<<(std::ostream & out, const AABB<n,T> & v) {
  out << '{' << v.min << ", " << v.max << "}";
  return out;
}

}
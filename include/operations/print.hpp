#pragma once

#include <iostream>

template < uint32_t n, typename T >
std::ostream& operator<<(std::ostream& out, const vec<n,T> & v) {
  out << '{';
  for (uint32_t i = 0; i < n; i++) {
    out << v[i];
    if (i != n-1) out << ',';
  }
  out << '}';

  return out;
}

template < uint32_t m, uint32_t n, typename T >
std::ostream& operator<<(std::ostream& out, const mat<m,n,T> & A) {
  out << '{';
  for (uint32_t i = 0; i < m; i++) {
    out << '{';
    for (uint32_t j = 0; j < n; j++) {
      out << A[i][j];
      if (j != n-1) out << ',';
    }
    out << '}';
    if (i != n-1) out << ',';
  }
  out << '}';

  return out;
}
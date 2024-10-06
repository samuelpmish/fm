#pragma once

#include "fm/types/matrix.hpp"

namespace fm {

template <Kind kind, u32 rows, u32 cols, typename T>
__host__ __device__ constexpr auto transpose(const matrix<kind, rows, cols, T>& A) {
  if constexpr (kind == Kind::Isotropic || kind == Kind::Diagonal ||
                kind == Kind::Symmetric) {
    return A;
  }

  if constexpr (kind == Kind::General) {
    matrix<kind, cols, rows, T> AT{};
    for (int i = 0; i < cols; i++) {
      for (int j = 0; j < rows; j++) {
        AT(i, j) = A(j, i);
      }
    }
    return AT;
  }

  if constexpr (kind == Kind::Skew) {
    return -A;
  }

  if constexpr (kind == Kind::Rotation) {
    return matrix<kind, rows, cols, T>{A.c, -A.s};
  }
}

}  // namespace fm
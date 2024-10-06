#pragma once

namespace fm {

/// frobenius norm of a matrix
template < Kind kind, u32 dim, typename T >
__host__ __device__ constexpr auto norm(const matrix<kind, dim, dim, T> & A) {

  if constexpr (kind == Kind::General || kind == Kind::Symmetric) {
    T total{};
    for (u32 i = 0; i < dim; i++) {
      for (u32 j = 0; j < dim; j++) {
        total += A(i,j) * A(i,j);
      }
    }
    return sqrt(total);
  }

  if constexpr (kind == Kind::Diagonal) {
    T total{};
    for (u32 i = 0; i < dim; i++) {
      total += A.data[i] * A.data[i];
    }
    return sqrt(total);
  }

  if constexpr (kind == Kind::Rotation) {
    if constexpr (dim == 2) return 1.414213562373095;
    if constexpr (dim == 3) {
      static_assert(always_false<T>{}, "unimplemented");
    }
  }

  if constexpr (kind == Kind::Skew) {
    if constexpr (dim == 2) return 1.414213562373095 * A.data;
    if constexpr (dim == 3) return 1.414213562373095 * norm(A.data);
  }

  if constexpr (kind == Kind::Isotropic) {
    return sqrt(dim * (A.data * A.data));
  }

}

}
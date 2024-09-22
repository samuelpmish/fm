#pragma once

#include "fm/operations/transpose.hpp"
#include "fm/operations/invariants.hpp"

namespace fm {

template < Kind kind, u32 dim, typename T >
constexpr auto inv(const matrix<kind, dim, dim, T> & A) {

  using matrix_type = matrix<kind, dim, dim, T>;

  if constexpr (kind == Kind::General) {

    if constexpr (dim == 1) {
      return mat1{1.0 / A(0, 0)};
    }

    if constexpr (dim == 2) {
      T inv_detA(1.0 / det(A));

      return matrix_type{{
        { A[1][1] * inv_detA, -A[0][1] * inv_detA},
        {-A[1][0] * inv_detA,  A[0][0] * inv_detA}
      }};
    }

    if constexpr (dim == 3) {
      auto inv_detA = 1.0 / det(A);

      matrix_type invA{};
      invA[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) * inv_detA;
      invA[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) * inv_detA;
      invA[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * inv_detA;
      invA[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) * inv_detA;
      invA[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) * inv_detA;
      invA[1][2] = (A[0][2] * A[1][0] - A[0][0] * A[1][2]) * inv_detA;
      invA[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]) * inv_detA;
      invA[2][1] = (A[0][1] * A[2][0] - A[0][0] * A[2][1]) * inv_detA;
      invA[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) * inv_detA;
      return invA;   
    }

  }

  if constexpr (kind == Kind::Symmetric) {

    if constexpr (dim == 2) {
      T inv_detA(1.0 / det(A));

      matrix_type invA;
      invA(0,0) =  A(1,1) * inv_detA;
      invA(1,0) = -A(0,1) * inv_detA;
      invA(1,1) =  A(0,0) * inv_detA;
      return invA;
    }

    if constexpr (dim == 3) {
      auto inv_detA = 1.0 / det(A);

      matrix_type invA{};
      invA(0,0) = (A(1,1) * A(2,2) - A(1,2) * A(2,1)) * inv_detA;
      invA(0,1) = (A(0,2) * A(2,1) - A(0,1) * A(2,2)) * inv_detA;
      invA(0,2) = (A(0,1) * A(1,2) - A(0,2) * A(1,1)) * inv_detA;
      invA(1,1) = (A(0,0) * A(2,2) - A(0,2) * A(2,0)) * inv_detA;
      invA(1,2) = (A(0,2) * A(1,0) - A(0,0) * A(1,2)) * inv_detA;
      invA(2,2) = (A(0,0) * A(1,1) - A(0,1) * A(1,0)) * inv_detA;
      return invA;   
    }

  }

  if constexpr (kind == Kind::Diagonal) {
    matrix_type invA{};
    for (u32 i = 0; i < dim; i++) {
      invA.data[i] = 1.0 / A.data[i];
    }
    return invA;
  }

  if constexpr (kind == Kind::Skew) {
    if constexpr (dim == 2) {
      return matrix_type{-1.0 / A.data};
    }
    if constexpr (dim == 3) {
      static_assert(always_false<T>{}, "Invalid operation: cannot invert a singular matrix");
    }
  }

  if constexpr (kind == Kind::Rotation) {
    return transpose(A);
  }

  if constexpr (kind == Kind::Isotropic) {
    return matrix_type{1.0 / A.data};
  }

}

}
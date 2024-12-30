#pragma once

#include "fm/operations/transpose.hpp"

namespace fm {
  
constexpr auto adj(const double & A) { return 1.0; }

template < Kind kind, u32 dim, typename T >
__host__ __device__ constexpr auto adj(const matrix<kind, dim, dim, T> & A) {

  using matrix_type = matrix<kind, dim, dim, T>;

  if constexpr (kind == Kind::General) {

    if constexpr (dim == 1) {
      return mat1{1.0};
    }

    if constexpr (dim == 2) {
      return matrix_type{{
        { A[1][1], -A[0][1]},
        {-A[1][0],  A[0][0]}
      }};
    }

    if constexpr (dim == 3) {
      matrix_type adjA{};
      adjA[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]);
      adjA[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]);
      adjA[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]);
      adjA[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]);
      adjA[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]);
      adjA[1][2] = (A[0][2] * A[1][0] - A[0][0] * A[1][2]);
      adjA[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
      adjA[2][1] = (A[0][1] * A[2][0] - A[0][0] * A[2][1]);
      adjA[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]);
      return adjA;   
    }

  }

  if constexpr (kind == Kind::Symmetric) {

    if constexpr (dim == 2) {
      matrix_type adjA;
      adjA(0,0) =  A(1,1);
      adjA(1,0) = -A(0,1);
      adjA(1,1) =  A(0,0);
      return adjA;
    }

    if constexpr (dim == 3) {
      matrix_type adjA{};
      adjA(0,0) = (A(1,1) * A(2,2) - A(1,2) * A(2,1));
      adjA(0,1) = (A(0,2) * A(2,1) - A(0,1) * A(2,2));
      adjA(0,2) = (A(0,1) * A(1,2) - A(0,2) * A(1,1));
      adjA(1,1) = (A(0,0) * A(2,2) - A(0,2) * A(2,0));
      adjA(1,2) = (A(0,2) * A(1,0) - A(0,0) * A(1,2));
      adjA(2,2) = (A(0,0) * A(1,1) - A(0,1) * A(1,0));
      return adjA;   
    }

  }

  if constexpr (kind == Kind::Diagonal) {
    auto detA = det(A);
    matrix_type adjA{};
    for (u32 i = 0; i < dim; i++) {
      adjA.data[i] = detA / A.data[i];
    }
    return adjA;
  }

  if constexpr (kind == Kind::Skew) {
    if constexpr (dim == 2) {
      return matrix_type{-A.data};
    }
    if constexpr (dim == 3) {
      return outer(A.data, A.data);
    }
  }

  if constexpr (kind == Kind::Rotation) {
    return transpose(A);
  }

  if constexpr (kind == Kind::Isotropic) {
    return matrix_type{det(A) / A.data};
  }

}

}
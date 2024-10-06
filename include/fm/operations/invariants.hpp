#pragma once

namespace fm {

__host__ __device__ constexpr auto tr(const double & A) {
  return A;
}

template < Kind kind, u32 dim, typename T >
__host__ __device__ constexpr auto tr(const matrix<kind, dim, dim, T> & A) {

  if constexpr (kind == Kind::General || kind == Kind::Symmetric) {
    if constexpr (dim == 1) return A(0,0);
    if constexpr (dim == 2) return A(0,0) + A(1,1);
    if constexpr (dim == 3) return A(0,0) + A(1,1) + A(2,2);
  }

  if constexpr (kind == Kind::Diagonal) {
    if constexpr (dim == 1) return A.data[0];
    if constexpr (dim == 2) return A.data[0] + A.data[1];
    if constexpr (dim == 3) return A.data[0] + A.data[1] + A.data[2];
  }

  if constexpr (kind == Kind::Rotation) {
    if constexpr (dim == 2) return 2 * A.c;
    if constexpr (dim == 3) return 3 - 4 * dot(A.s, A.s);
  }

  if constexpr (kind == Kind::Skew) {
    return T{};
  }

  if constexpr (kind == Kind::Isotropic) {
    return dim * A.data;
  }

}

__host__ __device__ constexpr auto det(const double & A) {
  return A;
}

template < Kind kind, u32 dim, typename T >
__host__ __device__ constexpr auto det(const matrix<kind, dim, dim, T> & A) {

  if constexpr (kind == Kind::General || kind == Kind::Symmetric) {
    if constexpr (dim == 1) {
      return A(0,0);
    }
    if constexpr (dim == 2) {
      return A(0,0) * A(1,1) - A(0,1) * A(1,0);
    }
    if constexpr (dim == 3) {
      return A(0,0) * A(1,1) * A(2,2) + A(0,1) * A(1,2) * A(2,0) +
             A(0,2) * A(1,0) * A(2,1) - A(0,0) * A(1,2) * A(2,1) -
             A(0,1) * A(1,0) * A(2,2) - A(0,2) * A(1,1) * A(2,0);
    }
  }

  if constexpr (kind == Kind::Diagonal) {
    if constexpr (dim == 2) return A.data[0] * A.data[1];
    if constexpr (dim == 3) return A.data[0] * A.data[1] * A.data[2];
  }

  if constexpr (kind == Kind::Isotropic) {
    if constexpr (dim == 2) return A.data * A.data;
    if constexpr (dim == 3) return A.data * A.data * A.data;
  }

  if constexpr (kind == Kind::Skew) { 
    if constexpr (dim == 2) return A.data * A.data;
    if constexpr (dim == 3) return T(0);
  }

  if constexpr (kind == Kind::Rotation) { return T(1); }

}

template < u32 which, Kind kind, u32 dim, typename T >
__host__ __device__ constexpr auto invariant(const matrix<kind, dim, dim, T> & A) {

  static_assert(dim == 2 || dim == 3, "invariants only defined for 2x2 and 3x3 matrices");

  if constexpr (which == 1) { return tr(A); }

  if constexpr (which == dim) { return det(A); }

  if constexpr (which == 2 && dim == 3) { 

    if constexpr (kind == Kind::General) {
      return A(0, 1) * A(1, 0) - A(0, 0) * A(1, 1) + 
             A(0, 2) * A(2, 0) + A(1, 2) * A(2, 1) - 
             A(0, 0) * A(2, 2) - A(1, 1) * A(2, 2);
    }

    if constexpr (kind == Kind::Symmetric) {
      return A.data[1] * A.data[1] + A.data[2] * A.data[2] - 
             A.data[0] * A.data[3] + A.data[4] * A.data[4] - 
             A.data[0] * A.data[5] - A.data[3] * A.data[5];
    }

    if constexpr (kind == Kind::Rotation) {
      return -tr(A);
    }

    if constexpr (kind == Kind::Skew) {
      return -dot(A.data, A.data);
    }

    if constexpr (kind == Kind::Diagonal) {
      return -A.data[0] * A.data[1] - A.data[0] * A.data[2] - A.data[1] * A.data[2];
    }

    if constexpr (kind == Kind::Isotropic) {
      return -3 * A.data * A.data;
    }

  }

}

}
#pragma once

namespace fm {

template < Kind kind, u32 dim, typename T >
constexpr auto tr(const matrix<kind, dim, dim, T> & A) {

  if constexpr (kind == Kind::General || kind == Kind::Symmetric) {
    if constexpr (dim == 2) return A(0,0) + A(1,1);
    if constexpr (dim == 3) return A(0,0) + A(1,1) + A(2,2);
  }

  if constexpr (kind == Kind::Diagonal) {
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

template < Kind kind, u32 dim, typename T >
constexpr auto det(const matrix<kind, dim, dim, T> & A) {

  if constexpr (kind == Kind::General || kind == Kind::Symmetric) {
    if constexpr (dim == 2) {
      return A(0,0) * A(1,1) - A(0,1) * A(1,0);
    }
    if constexpr (dim == 3) {
      return A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] +
             A[0][2] * A[1][0] * A[2][1] - A[0][0] * A[1][2] * A[2][1] -
             A[0][1] * A[1][0] * A[2][2] - A[0][2] * A[1][1] * A[2][0];
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
    if constexpr (dim == 2) return -A.data * A.data;
    if constexpr (dim == 3) return T(0);
  }

  if constexpr (kind == Kind::Rotation) { return T(1); }

}

template < u32 which, Kind kind, u32 dim, typename T >
constexpr auto invariant(const matrix<kind, dim, dim, T> & A) {

  if constexpr (which == 1) { return tr(A); }

  if constexpr (which == dim) { return det(A); }

  if constexpr (which == 2 && dim == 3) { 

    if constexpr (kind == Kind::General) {
      static_assert(always_false<T>{}, "unimplemented");
    }

    if constexpr (kind == Kind::Symmetric) {
      static_assert(always_false<T>{}, "unimplemented");
    }

    if constexpr (kind == Kind::Rotation) {
      static_assert(always_false<T>{}, "unimplemented");
    }

    if constexpr (kind == Kind::Skew) {
      static_assert(always_false<T>{}, "unimplemented");
    }

    if constexpr (kind == Kind::Diagonal) {
      static_assert(always_false<T>{}, "unimplemented");
    }

    if constexpr (kind == Kind::Isotropic) {
      static_assert(always_false<T>{}, "unimplemented");
    }

  }

}

}
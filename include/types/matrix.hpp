#pragma once

#include "type_aliases.hpp"

#include "types/vec.hpp"

namespace fm {

template < typename T >
struct always_false : std::false_type {};

enum class Kind { Isotropic, Diagonal, Skew, Rotation, Symmetric, General };

// common interface:
// static constexpr member nrows()
// static constexpr member ncols()
// constexpr accessor A(i,j) for reading the (i,j) element
template <Kind k, u32 rows, u32 cols = rows, typename T = double >
struct matrix;

///////////////////////////////////////////////////////////////////////////////

template <u32 dim, typename T>
struct matrix<Kind::Isotropic, dim, dim, T> {
  using type = T;
  static constexpr u32 nrows() { return dim; }
  static constexpr u32 ncols() { return dim; }
  constexpr const auto operator()(u32 i, u32 j) const { return (i == j) * data; }
  T data;
};
template < u32 dim, typename T >
using iso = matrix<Kind::Isotropic, dim, dim, T>;
using iso2 = matrix<Kind::Isotropic, 2, 2, double>;
using iso3 = matrix<Kind::Isotropic, 3, 3, double>;

///////////////////////////////////////////////////////////////////////////////

template <u32 dim, typename T>
struct matrix<Kind::Diagonal, dim, dim, T> {
  using type = T;
  static constexpr u32 nrows() { return dim; }
  static constexpr u32 ncols() { return dim; }
  constexpr const auto operator()(u32 i, u32 j) const { return (i == j) * data[i]; }
  vec<dim,T> data;
};
template < u32 dim, typename T >
using diag = matrix<Kind::Diagonal, dim, dim, T>;
using diag2 = matrix<Kind::Diagonal, 2, 2, double>;
using diag3 = matrix<Kind::Diagonal, 3, 3, double>;

///////////////////////////////////////////////////////////////////////////////

template <typename T>
struct matrix<Kind::Skew, 2, 2, T> {
  using type = T;
  static constexpr u32 nrows() { return 2; }
  static constexpr u32 ncols() { return 2; }
  const auto operator()(u32 i, u32 j) const { 
    return data * (i != j) * ((j > i) ? 1 : -1);
  }

  constexpr operator matrix<Kind::General, 2, 2, T>() {
    return {{{0.0, data}, {-data, 0.0}}};
  }

  T data;
};

template <typename T>
struct matrix<Kind::Skew, 3, 3, T> {
  using type = T;
  static constexpr u32 nrows() { return 3; }
  static constexpr u32 ncols() { return 3; }
  constexpr const auto operator()(u32 i, u32 j) const { 
    return abs(i32(i) - i32(j)) * data[0 * (i == 0) & (j == 1) + 
                                       1 * (i == 1) & (j == 2) +
                                       2 * (i == 1) & (j == 2)];
  }
  T data[3];
};

template < u32 dim, typename T >
using skew = matrix<Kind::Skew, dim, dim, T>;
using skew2 = matrix<Kind::Skew, 2, 2, double>;
using skew3 = matrix<Kind::Skew, 3, 3, double>;

///////////////////////////////////////////////////////////////////////////////

template <typename T>
struct matrix<Kind::Rotation, 2, 2, T> {
  using type = T;
  static constexpr u32 nrows() { return 2; }
  static constexpr u32 ncols() { return 2; }
  constexpr const auto operator()(u32 i, u32 j) const { 
    return (i == j) * c + (i32(i) - i32(j)) * s;
  }

  constexpr operator matrix<Kind::General, 2, 2, T>() {
    return {{{c, -s}, {s, c}}};
  }

  T c;
  T s;
};

template <typename T>
struct matrix<Kind::Rotation, 3, 3, T> {
  using type = T;
  static constexpr u32 nrows() { return 3; }
  static constexpr u32 ncols() { return 3; }
  constexpr const auto & operator()(u32 i, u32 j) const { 
    static_assert(always_false<T>{}, "unimplemented");
  }

  constexpr operator matrix<Kind::General, 3, 3, T>() {
    return {{
      {1-2*(s[1]*s[1]+s[2]*s[2]),   2*(s[0]*s[1]-   c*s[2]),   2*(s[0]*s[2]+   c*s[1])},
      {  2*(s[0]*s[1]+   c*s[2]), 1-2*(s[0]*s[0]+s[2]*s[2]),   2*(s[1]*s[2]-   c*s[0])},
      {  2*(s[0]*s[2]-   c*s[1]),   2*(s[1]*s[2]+   c*s[0]), 1-2*(s[0]*s[0]+s[1]*s[1])}
    }};
  }

  T c;
  T s[3];
};

using rot2 = matrix<Kind::Rotation, 2, 2, double>;
using rot3 = matrix<Kind::Rotation, 3, 3, double>;

template < typename T >
matrix<Kind::Rotation, 2, 2, T> RotationMatrix(T angle) {
  return {cos(angle * 0.5), sin(angle * 0.5)};
}

///////////////////////////////////////////////////////////////////////////////

template <u32 dim, typename T>
struct matrix<Kind::Symmetric, dim, dim, T> {
  using type = T;
  static constexpr u32 nrows() { return dim; }
  static constexpr u32 ncols() { return dim; }
  static constexpr u32 num_values = (dim*(dim+1)) / 2;

  constexpr auto & operator()(u32 i, u32 j) { return data[index(i,j)]; }
  constexpr const auto & operator()(u32 i, u32 j) const { return data[index(i,j)]; }

  constexpr u32 index(u32 i, u32 j) const {
    u32 i_upper = std::min(i, j);
    u32 j_upper = std::max(i, j);
    return j_upper + ((2 * dim - i_upper - 1) * i_upper) / 2;
  }

  T data[num_values];
};

template < u32 n, typename T >
using sym = matrix<Kind::Symmetric, n, n, T>;
using sym2 = matrix<Kind::Symmetric, 2, 2, double>;
using sym3 = matrix<Kind::Symmetric, 3, 3, double>;

///////////////////////////////////////////////////////////////////////////////

template <u32 rows, u32 cols, typename T>
struct matrix<Kind::General, rows, cols, T> {
  using type = T;
  static constexpr u32 nrows() { return rows; }
  static constexpr u32 ncols() { return cols; }

  constexpr auto & operator()(u32 i, u32 j) { return data[i][j]; }
  constexpr const auto & operator()(u32 i, u32 j) const { return data[i][j]; }

  constexpr auto & operator[](u32 i) { return data[i]; }
  constexpr const auto & operator[](u32 i) const { return data[i]; }

  vec<cols, T> data[rows];
};

template < u32 m, u32 n, typename T >
using mat = matrix<Kind::General, m, n, T>;
using mat2 = matrix<Kind::General, 2, 2, double>;
using mat3 = matrix<Kind::General, 3, 3, double>;

}  // namespace fm

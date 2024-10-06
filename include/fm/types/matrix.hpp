#pragma once

#include "fm/type_aliases.hpp"

#include "fm/types/vec.hpp"
#include "fm/macros.hpp"

#include <iostream>
#include <algorithm>

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
  __host__ __device__ static constexpr u32 nrows() { return dim; }
  __host__ __device__ static constexpr u32 ncols() { return dim; }

  __host__ __device__ constexpr const auto operator[](u32 i) const { 
    vec<dim,T> output{};
    output[i] = data;
    return output;
  }

  __host__ __device__ constexpr const auto operator()(u32 i) const { 
    vec<dim,T> output{};
    output[i] = data;
    return output;
  }

  __host__ __device__ constexpr const auto operator()(u32 i, u32 j) const { return (i == j) * data; }

  __host__ __device__ constexpr operator matrix<Kind::General, dim, dim, T>() const {
    matrix<Kind::General, dim, dim, T> output{};
    for (u32 i = 0; i < dim; i++) {
      output(i,i) = data;
    }
    return output;
  }

  T data;
};
template < u32 dim, typename T=double >
using iso = matrix<Kind::Isotropic, dim, dim, T>;
using iso2 = matrix<Kind::Isotropic, 2, 2, double>;
using iso3 = matrix<Kind::Isotropic, 3, 3, double>;

///////////////////////////////////////////////////////////////////////////////

template <u32 dim, typename T>
struct matrix<Kind::Diagonal, dim, dim, T> {
  using type = T;
  __host__ __device__ static constexpr u32 nrows() { return dim; }
  __host__ __device__ static constexpr u32 ncols() { return dim; }

  __host__ __device__ constexpr const auto operator()(u32 i, u32 j) const { return (i == j) * data[i]; }

  __host__ __device__ constexpr operator matrix<Kind::General, dim, dim, T>() const {
    matrix<Kind::General, dim, dim, T> output{};
    for (u32 i = 0; i < dim; i++) {
      output(i,i) = data[i];
    }
    return output;
  }

  vec<dim,T> data;
};
template < u32 dim, typename T=double >
using diag = matrix<Kind::Diagonal, dim, dim, T>;
using diag2 = matrix<Kind::Diagonal, 2, 2, double>;
using diag3 = matrix<Kind::Diagonal, 3, 3, double>;

///////////////////////////////////////////////////////////////////////////////

template <typename T>
struct matrix<Kind::Skew, 2, 2, T> {
  using type = T;
  __host__ __device__ static constexpr u32 nrows() { return 2; }
  __host__ __device__ static constexpr u32 ncols() { return 2; }
  static constexpr u32 num_values = 1;

  __host__ __device__ constexpr const auto operator()(u32 i, u32 j) const { 
    return data * (i != j) * ((j < i) ? 1 : -1);
  }

  __host__ __device__ constexpr operator matrix<Kind::General, 2, 2, T>() const {
    return {{{0.0, -data}, {data, 0.0}}};
  }

  T data;
};

template <typename T>
struct matrix<Kind::Skew, 3, 3, T> {
  using type = T;
  __host__ __device__ static constexpr u32 nrows() { return 3; }
  __host__ __device__ static constexpr u32 ncols() { return 3; }
  static constexpr u32 num_values = 3;

  __host__ __device__ constexpr const auto operator()(u32 i, u32 j) const { 
    if (i == j) return 0.0;
    if (i == 0 && j == 1) return -data[2];
    if (i == 1 && j == 0) return  data[2];
    if (i == 2 && j == 0) return -data[1];
    if (i == 0 && j == 2) return  data[1];
    if (i == 1 && j == 2) return -data[0];
    if (i == 2 && j == 1) return  data[0];

    return T{};
  }

  __host__ __device__ constexpr operator matrix<Kind::General, 3, 3, T>() const {
    return {{
      {     0.0, -data[2], +data[1]}, 
      {+data[2],      0.0, -data[0]},
      {-data[1], +data[0],     0.0}
    }};
  }

  vec<3,T> data;
};

template < u32 dim, typename T=double >
using skew = matrix<Kind::Skew, dim, dim, T>;
using skew2 = matrix<Kind::Skew, 2, 2, double>;
using skew3 = matrix<Kind::Skew, 3, 3, double>;

///////////////////////////////////////////////////////////////////////////////

template <typename T>
struct matrix<Kind::Rotation, 2, 2, T> {
  using type = T;
  __host__ __device__ static constexpr u32 nrows() { return 2; }
  __host__ __device__ static constexpr u32 ncols() { return 2; }
  __host__ __device__ constexpr const auto operator()(u32 i, u32 j) const { 
    return (i == j) * c + (i32(i) - i32(j)) * s;
  }

  __host__ __device__ constexpr operator matrix<Kind::General, 2, 2, T>() const {
    return {{{c, -s}, {s, c}}};
  }

  T c;
  T s;
};

template <typename T>
struct matrix<Kind::Rotation, 3, 3, T> {
  using type = T;
  __host__ __device__ static constexpr u32 nrows() { return 3; }
  __host__ __device__ static constexpr u32 ncols() { return 3; }
  __host__ __device__ constexpr const auto operator()(u32 i, u32 j) const { 
    if (i == 0 && j == 0) return 1-2*(s[1]*s[1]+s[2]*s[2]);
    if (i == 0 && j == 1) return   2*(s[0]*s[1]-   c*s[2]);
    if (i == 0 && j == 2) return   2*(s[0]*s[2]+   c*s[1]);
    if (i == 1 && j == 0) return   2*(s[0]*s[1]+   c*s[2]);
    if (i == 1 && j == 1) return 1-2*(s[0]*s[0]+s[2]*s[2]);  
    if (i == 1 && j == 2) return   2*(s[1]*s[2]-   c*s[0]);
    if (i == 2 && j == 0) return   2*(s[0]*s[2]-   c*s[1]); 
    if (i == 2 && j == 1) return   2*(s[1]*s[2]+   c*s[0]);
    if (i == 2 && j == 2) return 1-2*(s[0]*s[0]+s[1]*s[1]);

    return T{};
  }

  __host__ __device__ constexpr operator matrix<Kind::General, 3, 3, T>() const {
    return {{
      {1-2*(s[1]*s[1]+s[2]*s[2]),   2*(s[0]*s[1]-   c*s[2]),   2*(s[0]*s[2]+   c*s[1])},
      {  2*(s[0]*s[1]+   c*s[2]), 1-2*(s[0]*s[0]+s[2]*s[2]),   2*(s[1]*s[2]-   c*s[0])},
      {  2*(s[0]*s[2]-   c*s[1]),   2*(s[1]*s[2]+   c*s[0]), 1-2*(s[0]*s[0]+s[1]*s[1])}
    }};
  }

  T c;
  vec<3,T> s;
};

template < u32 dim, typename T=double >
using rot = matrix<Kind::Rotation, dim, dim, T>;
using rot2 = matrix<Kind::Rotation, 2, 2, double>;
using rot3 = matrix<Kind::Rotation, 3, 3, double>;

template < typename T >
__host__ __device__ constexpr matrix<Kind::Rotation, 2, 2, T> RotationMatrix(T angle) {
  return {cos(angle), sin(angle)};
}

template < typename T >
__host__ __device__ constexpr rot<3,T> RotationMatrix(const vec<3,T> & axis_angle) {
  T theta = norm(axis_angle);
  vec<3,T> normalized_axis = axis_angle / theta;
  return {std::cos(theta / 2), std::sin(theta / 2) * normalized_axis};
}

///////////////////////////////////////////////////////////////////////////////

template <u32 dim, typename T>
struct matrix<Kind::Symmetric, dim, dim, T> {
  using type = T;
  __host__ __device__ static constexpr u32 nrows() { return dim; }
  __host__ __device__ static constexpr u32 ncols() { return dim; }

  static constexpr u32 num_values = (dim*(dim+1)) / 2;

  __host__ __device__ constexpr auto & operator()(u32 i, u32 j) { return data[index(i,j)]; }
  __host__ __device__ constexpr const auto & operator()(u32 i, u32 j) const { return data[index(i,j)]; }

  __host__ __device__ constexpr u32 index(u32 i, u32 j) const {
    #ifdef __CUDACC__
      u32 i_upper = (i < j) ? i : j;
      u32 j_upper = (i < j) ? j : i;
    #else
      u32 i_upper = std::min(i, j);
      u32 j_upper = std::max(i, j);
    #endif

    return j_upper + ((2 * dim - i_upper - 1) * i_upper) / 2;
  }

  __host__ __device__ constexpr operator matrix<Kind::General, dim, dim, T>() const {
    matrix<Kind::General, dim, dim, T> output;
    for (u32 i = 0; i < dim; i++) {
      for (u32 j = 0; j < dim; j++) {
        output(i,j) = data[index(i,j)];
      }
    }
    return output;
  }

  vec<num_values, T> data;
};

template < u32 n, typename T >
using sym = matrix<Kind::Symmetric, n, n, T>;
using sym2 = matrix<Kind::Symmetric, 2, 2, double>;
using sym3 = matrix<Kind::Symmetric, 3, 3, double>;

///////////////////////////////////////////////////////////////////////////////

template <u32 rows, u32 cols, typename T>
struct matrix<Kind::General, rows, cols, T> {
  using type = T;
  __host__ __device__ static constexpr u32 nrows() { return rows; }
  __host__ __device__ static constexpr u32 ncols() { return cols; }

  __host__ __device__ constexpr auto & operator()(u32 i, u32 j) { return data[i][j]; }
  __host__ __device__ constexpr const auto & operator()(u32 i, u32 j) const { return data[i][j]; }

  __host__ __device__ constexpr auto & operator[](u32 i) { return data[i]; }
  __host__ __device__ constexpr const auto & operator[](u32 i) const { return data[i]; }

  vec<cols, T> data[rows];
};

template < u32 m, u32 n, typename T = double >
using mat = matrix<Kind::General, m, n, T>;
using mat1 = matrix<Kind::General, 1, 1, double>;
using mat2 = matrix<Kind::General, 2, 2, double>;
using mat3 = matrix<Kind::General, 3, 3, double>;
using mat4 = matrix<Kind::General, 4, 4, double>;

using mat2x3 = matrix<Kind::General, 2, 3, double>;
using mat2x4 = matrix<Kind::General, 2, 4, double>;
using mat3x2 = matrix<Kind::General, 3, 2, double>;
using mat3x4 = matrix<Kind::General, 3, 4, double>;
using mat4x2 = matrix<Kind::General, 4, 2, double>;
using mat4x3 = matrix<Kind::General, 4, 3, double>;

using mat2f   = matrix<Kind::General, 2, 2, float>;
using mat3f   = matrix<Kind::General, 3, 3, float>;
using mat4f   = matrix<Kind::General, 4, 4, float>;

using mat2x3f = matrix<Kind::General, 2, 3, float>;
using mat2x4f = matrix<Kind::General, 2, 4, float>;
using mat3x2f = matrix<Kind::General, 3, 2, float>;
using mat3x4f = matrix<Kind::General, 3, 4, float>;
using mat4x2f = matrix<Kind::General, 4, 2, float>;
using mat4x3f = matrix<Kind::General, 4, 3, float>;

template < Kind k, uint32_t m, uint32_t n, typename T >
__host__ __device__ constexpr int tensor_rank(matrix<k,m,n,T>) { return 2; }

template <int dim, typename T=double>
__host__ __device__ constexpr iso<dim, T> Identity() {
  return iso<dim,T>{1.0};
}

template < typename T >
__host__ __device__ constexpr mat<3,3,T> to_3x3(const mat<2,2,T> & A) {
  return mat<3,3,T>{{
    {A[0][0], A[0][1], 0.0},
    {A[1][0], A[1][1], 0.0},
    {    0.0,     0.0, 0.0}
  }};
}
 
template < typename T >
__host__ __device__ constexpr mat<2,2,T> to_2x2(const mat<3,3,T> & A) {
  return mat<2,2,T>{{
    {A[0][0], A[0][1]},
    {A[1][0], A[1][1]}
  }};
}

template < typename T>
__host__ __device__ vec<3,T> cross(const mat<3,2,T> & A) { 
  return vec<3, T>{
    A(1,0)*A(2,1)-A(2,0)*A(1,1),
    A(2,0)*A(0,1)-A(0,0)*A(2,1),
    A(0,0)*A(1,1)-A(1,0)*A(0,1)
  };
}

template < typename T>
__host__ __device__ vec<2,T> cross(const mat<2,1,T> & A) { 
  return vec<2, T>{A(1,0), -A(0,0)};
}

template < Kind k, u32 m, u32 n, typename T >
std::ostream & operator<<(std::ostream & out, const matrix<k, m, n, T> & A) {
  out << '{' << '\n';
  for (int i = 0; i < m; i++) {
    out << "  {" << A(i,0);
    for (int j = 1; j < n; j++) {
      out << ", " << A(i,j);
    }
    out << '}' << '\n';
  }
  out << '}';
  return out;
}

template <typename S, typename T, u32 m, u32 n>
__host__ __device__ constexpr auto outer(const vec<m,S> & u, const vec<n,T> & v) {
  mat<m,n,decltype(S{} * T{})> uvT{};
  for (u32 i = 0; i < m; i++) {
    for (u32 j = 0; j < n; j++) {
      uvT[i][j] = u[i] * v[j];
    }
  }
  return uvT;
}

////////////////////////////////////////////////////////////////////////////////

template <Kind kind, u32 m, u32 n, typename T>
__host__ __device__ constexpr auto chain_rule(const matrix<kind, m, n, T> & df_dx, const T & dx) {

  using output_t = decltype(inner(T{}, T{}));

  if constexpr (kind == Kind::Isotropic) {
    return matrix<kind, m, n, output_t>{inner(df_dx.data, dx)};
  }

  if constexpr (kind == Kind::Diagonal) {
    matrix<kind, m, n, output_t> df{};
    for (u32 k = 0; k < n; k++) {
      df.data[k] = inner(df_dx.data[k], dx);
    }
    return df;
  }

  if constexpr (kind == Kind::Skew || kind == Kind::Symmetric) {
    matrix<kind, m, n, output_t> df{};
    for (u32 k = 0; k < decltype(df_dx)::num_values; k++) {
      df.data[k] = inner(df_dx.data[k], dx);
    }
    return df;
  }

  if constexpr (kind == Kind::Rotation) {
    matrix<kind, m, n, output_t> df{};
    if constexpr (m == 2) {
      df.s = inner(df_dx.s, dx);
      df.c = inner(df_dx.c, dx);
    }
    if constexpr (m == 3) {
      df.s[0] = inner(df_dx.s[0], dx);
      df.s[1] = inner(df_dx.s[1], dx);
      df.s[2] = inner(df_dx.s[2], dx);
      df.c = inner(df_dx.c, dx);
    }
    return df;
  }

  if constexpr (kind == Kind::General) {
    matrix<kind, m, n, output_t> df{};
    for (u32 i = 0; i < m; i++) {
      for (u32 j = 0; j < n; j++) {
        df(i,j) = inner(df_dx(i,j), dx);
      }
    }
    return df;
  }

}

}  // namespace fm

#include "fm/operations/dot.hpp"
#include "fm/operations/eig.hpp"
#include "fm/operations/invariants.hpp"
#include "fm/operations/inverse.hpp"
#include "fm/operations/linear_solve.hpp"
#include "fm/operations/norm.hpp"
#include "fm/operations/operator_overloads.hpp"
#include "fm/operations/transpose.hpp"

#pragma once

#include <cmath>
#include <iostream>
#include <algorithm> // for std::{min,max}

namespace fm {

template < uint32_t dim, typename T = double >
struct vec { 
  using data_type = T;
  static constexpr int dimension = dim;

  constexpr vec() : data{} {}

  constexpr vec(T a) {
    static_assert(dim == 1, "vec(T) only supported for vec<1,T>");
    data[0] = a;
  }

  constexpr vec(T a, T b) {
    static_assert(dim == 2, "vec(T, T) only supported for vec<2,T>");
    data[0] = a;
    data[1] = b;
  }

  constexpr vec(T a, T b, T c) {
    static_assert(dim == 3, "vec(T, T, T) only supported for vec<3,T>");
    data[0] = a;
    data[1] = b;
    data[2] = c;
  }

  constexpr vec(T a, T b, T c, T d) {
    static_assert(dim == 4, "vec(T, T, T, T) only supported for vec<4,T>");
    data[0] = a;
    data[1] = b;
    data[2] = c;
    data[3] = d;
  }

  // only allow implicit conversion to `T` when the vec is degenerate
  constexpr operator typename std::conditional<dim == 1, T, void>::type () const {
    return data[0];
  }

  constexpr T & operator[](uint32_t i) { return data[i]; }
  constexpr const T & operator[](uint32_t i) const { return data[i]; }

  constexpr T & operator()(uint32_t i) { return data[i]; }
  constexpr const T & operator()(uint32_t i) const { return data[i]; }

  T data[dim];

};

template < uint32_t dim, typename T = float >
vec(const T (&)[dim]) -> vec<dim, T>;

template < uint32_t dim >
using vecf = vec<dim, float>;

using vec1f = vec<1, float>;
using vec1 = vec<1, double>;

using vec2f = vec<2, float>;
using vec2 = vec<2, double>;

using vec3f = vec<3, float>;
using vec3 = vec<3, double>;

using vec4f = vec<4, float>;
using vec4 = vec<4, double>;

constexpr uint32_t dimension(double) { return 1; }

template < uint32_t dim, typename T >
constexpr uint32_t dimension(vec < dim, T >) { return dim; }

constexpr auto outer(const double & u, const double & v) {
  return u * v;
}

template <typename T, int n>
constexpr auto outer(const double & u, const vec<n,T> & v) {
  return u * v;
}

template <typename T, int n>
constexpr auto outer(const vec<n,T> & u, const double & v) {
  return u * v;
}

////////////////////////////////////////////////////////////////////////////////

template < uint32_t dim, typename S, typename T >
constexpr auto operator!=(const vec< dim, S > & u, const vec< dim, T > & v) {
  for (int i = 0; i < dim; i++) {
    if (u[i] != v[i]) return true;
  }
  return false;
}

template < uint32_t dim, typename S, typename T >
constexpr auto operator==(const vec< dim, S > & u, const vec< dim, T > & v) {
  for (uint32_t i = 0; i < dim; i++) {
    if (u[i] != v[i]) return false;
  }
  return true;
}

template < uint32_t n, typename T >
constexpr void operator*=(vec< n, T > & u, const T & v) {
  for (uint32_t i = 0; i < n; i++) u[i] *= v;
}

template < uint32_t n, typename T >
constexpr void operator+=(vec< n, T > & u, const vec< n, T > & v) {
  for (uint32_t i = 0; i < n; i++) u[i] += v[i];
}

template < uint32_t n, typename T >
constexpr void operator-=(vec< n, T > & u, const vec< n, T > & v) {
  for (uint32_t i = 0; i < n; i++) u[i] -= v[i];
}

template < uint32_t dim, typename S, typename T >
constexpr auto operator+(const vec< dim, S > & u, const vec< dim, T > & v) {
  vec< dim, decltype(S{} + T{}) > out{};
  for (uint32_t i = 0; i < dim; i++) {
    out[i] = u[i] + v[i];
  }
  return out; 
}

template < uint32_t dim, typename T >
constexpr auto operator-(const vec< dim, T > & v) {
  vec< dim, T > out{};
  for (uint32_t i = 0; i < dim; i++) {
    out[i] = -v[i];
  }
  return out; 
}

template < uint32_t dim, typename S, typename T >
constexpr auto operator-(const vec< dim, S > & u, const vec< dim, T > & v) {
  vec< dim, decltype(S{} - T{}) > out{};
  for (uint32_t i = 0; i < dim; i++) {
    out[i] = u[i] - v[i];
  }
  return out; 
}

template < uint32_t dim, typename S, typename T >
constexpr auto operator*(const vec< dim, S > & u, const vec< dim, T > & v) {
  vec< dim, decltype(S{} * T{}) > out{};
  for (uint32_t i = 0; i < dim; i++) {
    out[i] = u[i] * v[i];
  }
  return out; 
}

template < uint32_t dim, typename S, typename T >
constexpr auto operator/(const vec< dim, S > & u, const vec< dim, T > & v) {
  vec< dim, decltype(S{} / T{}) > out{};
  for (uint32_t i = 0; i < dim; i++) {
    out[i] = u[i] / v[i];
  }
  return out; 
}

////////////////////////////////////////////////////////////////////////////////

template < uint32_t dim, typename T >
constexpr auto operator*(const double & u, const vec< dim, T > & v) {
  vec< dim, T > out{};
  for (uint32_t i = 0; i < dim; i++) {
    out[i] = u * v[i];
  }
  return out; 
}

constexpr double dot(const double & u, const double & v) {
  return u * v;
}

template < uint32_t dim, typename T >
constexpr auto dot(const T & u, const vec< dim, T > & v) {
  vec< dim, T > out{};
  for (uint32_t i = 0; i < dim; i++) {
    out[i] = u * v[i];
  }
  return out; 
}

template < uint32_t dim, typename T >
constexpr auto operator*(const vec< dim, T > & v, const double & u) {
  vec< dim, T > out{};
  for (uint32_t i = 0; i < dim; i++) {
    out[i] = v[i] * u;
  }
  return out; 
}

template < uint32_t dim, typename T >
constexpr auto dot(const vec< dim, T > & v, const T & u) {
  vec< dim, T > out{};
  for (uint32_t i = 0; i < dim; i++) {
    out[i] = v[i] * u;
  }
  return out; 
}

template < uint32_t dim, typename T >
constexpr auto operator/(const double & u, const vec< dim, T > & v) {
  vec< dim, T > out{};
  for (uint32_t i = 0; i < dim; i++) {
    out[i] = u / v[i];
  }
  return out;
}

template < uint32_t dim, typename T >
constexpr auto operator/(const vec< dim, T > & v, const double & u) {
  vec< dim, T > out{};
  for (uint32_t i = 0; i < dim; i++) {
    out[i] = v[i] / u;
  }
  return out;
}

template < uint32_t dim, typename T >
constexpr auto norm_squared(const vec< dim, T > & v) {
  T out{};
  for (uint32_t i = 0; i < dim; i++) {
    out += v[i] * v[i];
  }
  return out; 
}

template < uint32_t dim, typename T >
constexpr auto norm(const vec< dim, T > & v) {
  T out{};
  for (uint32_t i = 0; i < dim; i++) {
    out += v[i] * v[i];
  }
  return std::sqrt(out); 
}

template < uint32_t dim, typename T >
constexpr auto normalize(const vec< dim, T > & v) {
  T norm_v = norm(v);
  if (norm_v != 0) {
    return v / norm_v;
  } else {
    return vec< dim, T >{};
  }
}

template < uint32_t dim, typename S, typename T >
constexpr auto dot(const vec< dim, S > & u, const vec< dim, T > & v) {
  decltype(S{} / T{}) total{};
  for (uint32_t i = 0; i < dim; i++) {
    total += u[i] * v[i];
  }
  return total;
}

// vector-vector products
template <typename T>
constexpr vec<3,T> cross(const vec<3,T> & u, const vec<3,T> & v) {
  return vec<3,T>{
    u[1] * v[2] - u[2] * v[1], 
    u[2] * v[0] - u[0] * v[2],
    u[0] * v[1] - u[1] * v[0] 
  };
}

////////////////////////////////////////////////////////////////////////////////

template < uint32_t n, typename T>
constexpr vec<n,T> abs(const vec<n,T> & v) {
  vec<n,T> absv{};
  for (int i = 0; i < n; i++) {
    absv[i] = std::abs(v[i]);
  }
  return absv;
}

template < uint32_t n, typename T >
constexpr T min(const vec<n,T> & v) { 
  T minval = v[0];
  for (int i = 0; i < n; i++) {
    minval = std::min(minval, v[i]);
  }
  return minval;
}

template < uint32_t n, typename T >
inline vec<n,T> min(const vec<n,T> & v, T value) {
  vec<n,T> out{};
  for (int i = 0; i < n; i++) {
    out[i] = std::min(v[i], value);
  }
  return out;
}

template < uint32_t n, typename T >
inline vec<n,T> min(T value, const vec<n,T> & v) { return min(v, value); }

template < uint32_t n, typename T >
constexpr T max(const vec<n,T> & v) { 
  T maxval = v[0];
  for (int i = 0; i < n; i++) {
    maxval = std::max(maxval, v[i]);
  }
  return maxval;
}

template < uint32_t n, typename T >
constexpr vec<n,T> max(const vec<n,T> & v, T value) {
  vec<n,T> out{};
  for (int i = 0; i < n; i++) {
    out[i] = std::max(v[i], value);
  }
  return out;
}

template < int n, typename T >
constexpr vec<n,T> max(T value, const vec<n,T> & v) { return max(v, value); }

template < typename T >
T clamp(const T & value, const T & lower, const T & upper) {
  return std::max(lower, std::min(value, upper));
}

template < uint32_t n, typename T >
inline vec<n,T> clamp(const vec<n,T> & v, T lower, T upper) {
  vec<n,T> out{};
  for (uint32_t i = 0; i < n; i++) {
    out[i] = clamp(v[i], lower, upper);
  }
  return out;
}

////////////////////////////////////////////////////////////////////////////////

template < uint32_t n, typename T >
T inner(const vec<n, T> & u, const vec<n, T> & v) {
  return dot(u, v);
}

template < uint32_t n, typename T >
auto chain_rule(const vec<n, T> & df_dx, T dx) {
  vec<n, decltype(inner(T{}, T{}))> df{};
  for (int i = 0; i < n; i++) {
    df[i] = inner(df_dx[i], dx);
  }
  return df;
}

////////////////////////////////////////////////////////////////////////////////

template < typename T >
inline vec<2,T> xy(const vec<3,T> & v) { return {v[0], v[1]}; }

template < typename T >
inline vec<2,T> xy(const vec<2,T> & v) { return v; }

template < typename T >
inline vec<3,T> xyz(const vec<2,T> & v) { return {v[0], v[1], T{}}; }

template < typename T >
inline vec<3,T> xyz(const vec<3,T> & v) { return v; }

////////////////////////////////////////////////////////////////////////////////

constexpr int tensor_rank(double) { return 0; }

template < int n, typename T >
constexpr int tensor_rank(vec<n,T>) { return 1; }

template < int n >
std::ostream & operator<<(std::ostream & out, vec< n, double > v) {
  out << '{';
  for (int i = 0; i < n; i++) {
    out << v[i];
    if (i != n - 1) out << ", ";
  }
  out << '}';
  return out;
}

////////////////////////////////////////////////////////////////////////////////

template < typename T, uint32_t n >
std::ostream & operator<<(std::ostream & out, vec<n,T> v) {
  out << '{' << v(0);
  for (int i = 1; i < n; i++) {
    out << ", " << v(i);
  }
  out << '}';
  return out;
}

}
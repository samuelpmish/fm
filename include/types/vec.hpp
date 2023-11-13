#pragma once

#include <cmath>

template < uint32_t dim, typename T = double >
struct vec { 
  using data_type = T;
  static constexpr int dimension = dim;

  T data[dim];

  T & operator[](uint32_t i) { return data[i]; }
  const T & operator[](uint32_t i) const { return data[i]; }

  T & operator()(uint32_t i) { return data[i]; }
  const T & operator()(uint32_t i) const { return data[i]; }
};

template < uint32_t dim, typename T = float >
vec(const T (&)[dim]) -> vec<dim, T>;

using vec2f = vec<2, float>;
using vec2 = vec<2, double>;

using vec3f = vec<3, float>;
using vec3 = vec<3, double>;

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



////
template < uint32_t dim, typename T >
constexpr auto operator*(const double & u, const vec< dim, T > & v) {
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
constexpr auto norm(const vec< dim, T > & v) {
  T out{};
  for (uint32_t i = 0; i < dim; i++) {
    out += v[i] * v[i];
  }
  return std::sqrt(out); 
}

template < uint32_t dim, typename T >
constexpr auto normalize(const vec< dim, T > & v) {
  return (v / norm(v));
}

template < uint32_t dim, typename S, typename T >
constexpr auto dot(const vec< dim, S > & u, const vec< dim, T > & v) {
  decltype(S{} / T{}) total{};
  for (uint32_t i = 0; i < dim; i++) {
    total += u[i] * v[i];
  }
  return total;
}
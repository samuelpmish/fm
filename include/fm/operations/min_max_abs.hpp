#pragma once

#include "fm/types/vec.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace fm {

template < uint32_t n, typename T>
__host__ __device__ constexpr vec<n,T> abs(const vec<n,T> & v) {
  vec<n,T> absv{};
  for (int i = 0; i < n; i++) {
    absv[i] = std::abs(v[i]);
  }
  return absv;
}

__host__ __device__ constexpr float min(const float a, const float b) { 
  return (a < b) ? a : b;
}

__host__ __device__ constexpr double min(const double a, const double b) { 
  return (a < b) ? a : b;
}

__host__ __device__ constexpr float max(const float a, const float b) { 
  return (a < b) ? b : a;
}

__host__ __device__ constexpr double max(const double a, const double b) { 
  return (a < b) ? b : a;
}

template < uint32_t n, typename T >
__host__ __device__ constexpr T min(const vec<n,T> & v) { 
  T minval = v[0];
  for (int i = 0; i < n; i++) {
    minval = fm::min(minval, v[i]);
  }
  return minval;
}

template < uint32_t n, typename T >
__host__ __device__ constexpr vec<n,T> min(const vec<n,T> & v, T value) {
  vec<n,T> out{};
  for (int i = 0; i < n; i++) {
    out[i] = fm::min(v[i], value);
  }
  return out;
}

template < uint32_t n, typename T >
__host__ __device__ constexpr vec<n,T> min(T value, const vec<n,T> & v) { return min(v, value); }

template < uint32_t n, typename T >
__host__ __device__ constexpr vec<n,T> min(const vec<n,T> & u, const vec<n,T> & v) {
  vec<n,T> out{};
  for (int i = 0; i < n; i++) {
    out[i] = fm::min(u[i], v[i]);
  }
  return out;
}

template < uint32_t n, typename T >
__host__ __device__ constexpr T max(const vec<n,T> & v) { 
  T maxval = v[0];
  for (int i = 0; i < n; i++) {
    maxval = fm::max(maxval, v[i]);
  }
  return maxval;
}

template < uint32_t n, typename T >
__host__ __device__ constexpr vec<n,T> max(const vec<n,T> & v, T value) {
  vec<n,T> out{};
  for (int i = 0; i < n; i++) {
    out[i] = fm::max(v[i], value);
  }
  return out;
}

template < int n, typename T >
__host__ __device__ constexpr vec<n,T> max(T value, const vec<n,T> & v) { return max(v, value); }

template < uint32_t n, typename T >
__host__ __device__ constexpr vec<n,T> max(const vec<n,T> & u, const vec<n,T> & v) {
  vec<n,T> out{};
  for (int i = 0; i < n; i++) {
    out[i] = fm::max(u[i], v[i]);
  }
  return out;
}

template < typename T >
__host__ __device__ constexpr T clamp(const T & value, const T & lower, const T & upper) {
  return fm::max(lower, fm::min(value, upper));
}

template < uint32_t n, typename T >
__host__ __device__ constexpr vec<n,T> clamp(const vec<n,T> & v, T lower, T upper) {
  vec<n,T> out{};
  for (uint32_t i = 0; i < n; i++) {
    out[i] = clamp(v[i], lower, upper);
  }
  return out;
}

}
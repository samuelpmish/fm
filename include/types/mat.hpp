#pragma once

#include <iostream>

template < uint32_t nrows, uint32_t ncols, typename T = double >
struct mat { 
  using data_type = T;
  static constexpr int dimensions[2] = {nrows,ncols};
  static constexpr fm::type type = fm::type::mat;

  vec<ncols, T> data[nrows];

  auto & operator[](uint32_t i) { return data[i]; }
  const auto & operator[](uint32_t i) const { return data[i]; }

  auto & operator()(uint32_t i) { return data[i]; }
  const auto & operator()(uint32_t i) const { return data[i]; }

  auto & operator()(uint32_t i, uint32_t j) { return data[i][j]; }
  const auto & operator()(uint32_t i, uint32_t j) const { return data[i][j]; }

//  operator sym<nrows,T>() {
//    sym<n,T> output;
//    for (uint32_t i = 0; i < n; i++) {
//      for (uint32_t j = i; j < n; j++) {
//        output(i,j) = (data[i][j] + data[j][i]) * 0.5;
//      } 
//    } 
//    return output;
//  }
};

using mat2f = mat<2, 2, float>;
using mat2 = mat<2, 2, double>;

using mat3f = mat<3, 3, float>;
using mat3 = mat<3, 3, double>;

using mat6 = mat<6, 6, double>;

template < uint32_t m, uint32_t n, typename S, typename T >
constexpr auto operator!=(const mat< m, n, S > & u, const mat< m, n, T > & v) {
  for (int i = 0; i < m; i++) {
    if (u[i] != v[i]) return true;
  }
  return false;
}

template < uint32_t m, uint32_t n, typename S, typename T >
constexpr auto operator==(const mat< m, n, S > & u, const mat< m, n, T > & v) {
  for (int i = 0; i < m; i++) {
    if (u[i] != v[i]) return false;
  }
  return true;
}

template <uint32_t m, uint32_t n, typename T>
constexpr auto transpose(const mat<m, n, T>& A) {
  mat<n, m, T> AT{};
  for (uint32_t i = 0; i < n; i++) {
    for (uint32_t j = 0; j < m; j++) {
      AT[i][j] = A[j][i];
    }
  }
  return AT;
}

template <typename T, uint32_t n>
constexpr T tr(const mat<n, n, T>& A) {
  T trA{};
  for (uint32_t i = 0; i < n; i++) {
    trA = trA + A(i,i);
  }
  return trA;
}
 
template <typename T, uint32_t n>
constexpr auto dev(const mat<n, n, T>& A) {
  auto devA = A;
  auto trA = tr(A);
  for (uint32_t i = 0; i < n; i++) {
    devA(i,i) -= trA / n;
  }
  return devA;
}

template < uint32_t m, uint32_t n, typename T >
constexpr void operator+=(mat< m, n, T > & A, const mat< m, n, T > & B) {
  for (uint32_t i = 0; i < m; i++) A[i] += B[i];
}

template < uint32_t m, uint32_t n, typename T >
constexpr void operator*=(mat< m, n, T > & A, const T & scale) {
  for (uint32_t i = 0; i < m; i++) A[i] *= scale;
}

template < uint32_t m, uint32_t n, typename S, typename T >
constexpr auto operator+(const mat< m, n, S > & A, const mat< m, n, T > & B) {
  mat< m, n, decltype(S{} + T{}) > out{};
  for (uint32_t i = 0; i < m; i++) {
    out[i] = A[i] + B[i];
  }
  return out; 
}

template < uint32_t m, uint32_t n, typename T >
constexpr auto operator*(const double & scale, const mat< m, n, T > & A) {
  mat< m, n, T > out{};
  for (uint32_t i = 0; i < m; i++) {
    out[i] = scale * A[i];
  }
  return out; 
}

template < uint32_t m, uint32_t n, typename T >
constexpr auto operator*(const mat< m, n, T > & A, const double & scale) {
  return scale * A;
}

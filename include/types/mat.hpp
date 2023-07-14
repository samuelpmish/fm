#pragma once

#include <iostream>

template < uint32_t nrows, uint32_t ncols, typename T = float >
struct mat { 
  vec<ncols, T> data[nrows];

  auto & operator[](uint32_t i) { return data[i]; }
  const auto & operator[](uint32_t i) const { return data[i]; }

  auto & operator()(uint32_t i) { return data[i]; }
  const auto & operator()(uint32_t i) const { return data[i]; }

  auto & operator()(uint32_t i, uint32_t j) { return data[i][j]; }
  const auto & operator()(uint32_t i, uint32_t j) const { return data[i][j]; }
};

using mat2f = mat<2, 2, float>;
using mat2 = mat<2, 2, double>;

using mat3f = mat<3, 3, float>;
using mat3 = mat<3, 3, double>;

using mat6 = mat<6, 6, double>;

template <uint32_t m, uint32_t n, typename T>
constexpr auto transpose(const mat<m, n, T>& A) {
  mat<n, m, T> AT{};
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      AT[i][j] = A[j][i];
    }
  }
  return AT;
}

template <typename T>
constexpr auto det(const mat<2, 2, T>& A) {
  return A[0][0] * A[1][1] - A[0][1] * A[1][0];
}

template <typename T>
constexpr auto det(const mat<3, 3, T>& A) {
  return A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] +
         A[0][2] * A[1][0] * A[2][1] - A[0][0] * A[1][2] * A[2][1] -
         A[0][1] * A[1][0] * A[2][2] - A[0][2] * A[1][1] * A[2][0];
}

template <typename T>
constexpr mat<2, 2, T> inv(const mat<2, 2, T>& A) {
  T inv_detA(1.0 / det(A));

  mat<2, 2, T> invA{};

  invA[0][0] =  A[1][1] * inv_detA;
  invA[0][1] = -A[0][1] * inv_detA;
  invA[1][0] = -A[1][0] * inv_detA;
  invA[1][1] =  A[0][0] * inv_detA;

  return invA;
}

template < typename T >
constexpr mat<3, 3, T> inv(const mat<3, 3, T>& A) {
  auto inv_detA = 1.0 / det(A);

  mat<3, 3, T> invA{};

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

template <typename T, int n>
constexpr T tr(const mat<n, n, T>& A) {
  T trA{};
  for (int i = 0; i < n; i++) {
    trA = trA + A(i,i);
  }
  return trA;
}
 
template <typename T, int n>
constexpr auto dev(const mat<n, n, T>& A) {
  auto devA = A;
  auto trA = tr(A);
  for (int i = 0; i < n; i++) {
    devA(i,i) -= trA / n;
  }
  return devA;
}

template < uint32_t m, uint32_t n, typename S, typename T >
constexpr auto operator+(const mat< m, n, S > & A, const mat< m, n, T > & B) {
  mat< m, n, decltype(S{} + T{}) > out{};
  for (int i = 0; i < m; i++) {
    out[i] = A[i] + B[i];
  }
  return out; 
}

template < uint32_t m, uint32_t n, typename T >
constexpr auto operator*(const double & scale, const mat< m, n, T > & A) {
  mat< m, n, T > out{};
  for (int i = 0; i < m; i++) {
    out[i] = scale * A[i];
  }
  return out; 
}

template < uint32_t m, uint32_t n, typename T >
constexpr auto operator*(const mat< m, n, T > & A, const double & scale) {
  return scale * A;
}

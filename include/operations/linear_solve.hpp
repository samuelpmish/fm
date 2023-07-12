#pragma once

#include "types/vec.hpp"
#include "types/mat.hpp"

template < uint32_t n, typename T >
constexpr vec<n,T> linear_solve(mat<n,n,T> A, vec<n,T> b) {

  constexpr auto abs = [](T x) { return (x < 0) ? -x : x; };

  vec<n,T> x{};

  for (int i = 0; i < n; i++) {
    // Search for maximum in this column
    T max_val = abs(A[i][i]);

    int max_row = i;
    for (int j = i + 1; j < n; j++) {
      if (abs(A(j, i)) > max_val) {
        max_val = abs(A[j][i]);
        max_row = j;
      }
    }

    // Swap maximum row with current row
    T tmp = b[max_row];
    b[max_row] = b[i];
    b[i] = tmp;
    for (int j = 0; j < n; j++) {
      tmp = A[max_row][j];
      A[max_row][j] = A[i][j];
      A[i][j] = tmp;
    }

    // zero entries below in this column
    for (int j = i + 1; j < n; j++) {
      T c = -A[j][i] / A[i][i];

      for (int k = i + 1; k < n; k++) {
        A[j][k] += c * A[i][k];
      }
      b[j] += c * b[i];
      A[j][i] = 0;
    }
  }

  // Solve equation Ax=b for an upper triangular matrix A
  for (int i = n - 1; i >= 0; i--) {
    x[i] = b[i] / A[i][i];
    for (int j = i - 1; j >= 0; j--) {
      b[j] -= A[j][i] * x[i];
    }
  }

  return x;
}

template <uint32_t n, typename T>
constexpr mat<n,n,T> inv(const mat<n,n,T>& A) {
  mat<n,n,T> invA{};
  for (int j = 0; j < n; j++) {
    vec<n,T> e_j{}; e_j[j] = 1.0;
    auto col = linear_solve(A, e_j);
    for (int i = 0; i < n; i++) {
      invA[i][j] = col[i];
    }
  }
  return invA;
}
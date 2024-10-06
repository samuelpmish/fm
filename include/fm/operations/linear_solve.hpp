#pragma once

#include "fm/types/vec.hpp"
#include "fm/types/matrix.hpp"

namespace fm {

template < Kind kind, uint32_t n, typename T >
__host__ __device__ constexpr auto linear_solve(const matrix<kind, n, n, T>& A, const vec<n,T>& b) {

  vec<n,T> x{};

  if constexpr (kind == Kind::General) {

    constexpr auto abs = [](T x) { return (x < 0) ? -x : x; };

    auto U = A;
    auto y = b;

    for (int i = 0; i < n; i++) {
      // Search for maximum in this column
      T max_val = abs(A[i][i]);

      int max_row = i;
      for (int j = i + 1; j < n; j++) {
        if (abs(U(j, i)) > max_val) {
          max_val = abs(U[j][i]);
          max_row = j;
        }
      }

      // Swap maximum row with current row
      T tmp = y[max_row];
      y[max_row] = y[i];
      y[i] = tmp;
      for (int j = 0; j < n; j++) {
        tmp = U[max_row][j];
        U[max_row][j] = U[i][j];
        U[i][j] = tmp;
      }

      // zero entries below in this column
      for (int j = i + 1; j < n; j++) {
        T c = -U[j][i] / U[i][i];

        for (int k = i + 1; k < n; k++) {
          U[j][k] += c * U[i][k];
        }
        y[j] += c * y[i];
        U[j][i] = 0;
      }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    for (int i = n - 1; i >= 0; i--) {
      x[i] = y[i] / U[i][i];
      for (int j = i - 1; j >= 0; j--) {
        y[j] -= U[j][i] * x[i];
      }
    }

  }

  if constexpr (kind == Kind::Symmetric) {

    auto U = A;

    // solve by LDLT factorization 
    for (int k = 0; k < n; k++) {
      x[k] = b[k];
      for (int i = k + 1; i < n; i++) {
        auto c = U(k,i)/U(k,k);
        for (int j = i; j < n; j++) {
          U(i, j) -= c * U(k, j);
        }
      }
    }

    for (int j = 0; j < n; j++) {
      auto s = 1.0 / U(j,j);
      x[j] *= s;
      for (int i = j+1; i < n; i++) {
        x[i] -= U(j,i) * x[j];
        U(i,j) *= s;
      }
      U(j,j) = 1.0;
    }

    for (int j = n - 1; j >= 0; j--) {
      for (int i = j - 1; i >= 0; i--) {
        x[i] -= U(i,j) * x[j];
      }
    }

  }

  if constexpr (kind == Kind::Diagonal) {
    for (int i = 0; i < n; i++) {
      x[i] = b[i] / A.data[i];
    } 
  }

  if constexpr (kind == Kind::Isotropic) {
    auto scale = 1.0 / A.data;
    for (int i = 0; i < n; i++) {
      x[i] = b[i] * scale;
    } 
  }

  if constexpr (kind == Kind::Rotation) {

  }

  return x;

}

}
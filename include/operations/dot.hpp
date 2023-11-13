#pragma once

#include <cinttypes>

template < typename S, typename T >
constexpr auto dot(const S & A, const T & B) {
  constexpr auto S_type = S::type;
  constexpr auto T_type = T::type;

  if constexpr (S_type == fm::type::vec) {

    // vec . vec
    if constexpr (T_type == fm::type::vec) {
      static_assert(S::dimension == T::dimension);
      constexpr int n = S::dimension;
      using return_type = decltype(A[0] * B[0]);
      return_type output{};
      for (int i = 0; i < n; i++) {
        output += A[i] * B[i];
      }
      return output;
    }

    // vec . {mat/sym}
    if constexpr (T_type == fm::type::mat || T_type == fm::type::sym) {
      static_assert(S::dimension == T::dimensions[0]);
      constexpr int m = S::dimension;
      constexpr int n = T::dimensions[1];
      using return_type = decltype(A[0] * B(0,0));
      vec< n, return_type > output{};
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
          output[i] += A[j] * B(j,i);
        }
      }
      return output;
    }

  }

  if constexpr (S_type == fm::type::mat || S_type == fm::type::sym) {

    // {mat/sym} . vec
    if constexpr (T_type == fm::type::vec) {
      static_assert(S::dimensions[1] == T::dimension);
      constexpr int m = S::dimensions[0];
      constexpr int n = S::dimensions[1];
      using return_type = decltype(A(0,0) * B[0]);
      vec< m, return_type > output{};
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
          output[i] += A(i,j) * B[j];
        }
      }
      return output;
    }

    // {mat/sym} . mat
    if constexpr (T_type == fm::type::mat) {
      static_assert(S::dimensions[1] == T::dimensions[0]);
      constexpr int m = S::dimensions[0];
      constexpr int n = T::dimensions[0];
      constexpr int p = T::dimensions[1];
      using return_type = decltype(A(0,0) * B[0][0]);
      mat< m, p, return_type > output{};
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
          for (int k = 0; k < n; k++) {
            output[i][j] += A(i,k) * B[k][j];
          }
        }
      }
      return output;
    }

    // {mat/sym} . sym
    if constexpr (T_type == fm::type::sym) {
      static_assert(S::dimensions[1] == T::dimensions[0], "incompatible matrix dimensions for dot product");
      constexpr int m = S::dimensions[0];
      constexpr int n = T::dimensions[0];
      constexpr int p = T::dimensions[1];
      using return_type = decltype(A(0,0) * B(0,0));
      mat< m, p, return_type > output{};
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < p; j++) {
          for (int k = 0; k < m; k++) {
            output[i][j] += A(i,k) * B(k,j);
          }
        }
      }
      return output;
    }

  }

  //static_assert(false, "unsupported types");

}

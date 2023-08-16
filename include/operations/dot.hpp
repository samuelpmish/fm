#pragma once

#include <cinttypes>

#if 0
template < uint32_t dim, typename S, typename T >
constexpr auto dot(const vec< dim, S > & u, const vec< dim, T > & v) {
  decltype(S{} / T{}) out{};
  for (int i = 0; i < dim; i++) {
    out += u[i] * v[i];
  }
  return out; 
}

template <uint32_t m, uint32_t n, typename S, typename T>
constexpr auto dot(const vec<m,S>& A, const mat<m,n,T>& B)
{
  vec<n,decltype(S{} * T{})> AB{};
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      AB[i] = AB[i] + A[j] * B[j][i];
    }
  }
  return AB;
}

template <uint32_t m, uint32_t n, typename S, typename T>
constexpr auto dot(const mat<m,n,S>& A, const vec<n,T>& B) {
  vec<m, decltype(S{} * T{})> AB{};
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      AB[i] = AB[i] + A[i][j] * B[j];
    }
  }
  return AB;
}

template <uint32_t m, uint32_t n, uint32_t p, typename S, typename T>
constexpr auto dot(const mat<m,n,S>& A, const mat<n,p,T>& B) {
  mat<m, p, decltype(S{} * T{})> AB{};
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < p; j++) {
      for (int k = 0; k < n; k++) {
        AB[i][j] = AB[i][j] + A[i][k] * B[k][j];
      }
    }
  }
  return AB;
}

template <uint32_t n, typename S, typename T>
constexpr auto dot(const vec<n,S>& y, const sym<n,T>& A) {
  vec<n, decltype(S{} * T{})>  yA{};
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      yA[j] += y[i] * A(i,j);
    }
  }
  return yA;
}

template <uint32_t n, typename S, typename T>
constexpr auto dot(const sym<n,S>& A, const vec<n,T>& x) {
  vec<n, decltype(S{} * T{})>  Ax{};
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Ax[i] += A(i,j) * x[j];
    }
  }
  return Ax;
}

#else

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
#endif

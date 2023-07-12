#pragma once

template < uint32_t dim, typename S, typename T >
constexpr auto dot(const vec< dim, S > & u, const vec< dim, T > & v) {
  decltype(S{} / T{}) out{};
  for (int i = 0; i < dim; i++) {
    out += u[i] * v[i];
  }
  return out; 
}

template <uint32_t m, uint32_t n,typename S, typename T>
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

template <uint32_t m, uint32_t n,typename S, typename T>
constexpr auto dot(const mat<m,n,S>& A, const vec<n,T>& B)
{
  vec<m, decltype(S{} * T{})> AB{};
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      AB[i] = AB[i] + A[i][j] * B[j];
    }
  }
  return AB;
}

template <uint32_t m, uint32_t n, uint32_t p, typename S, typename T>
constexpr auto dot(const mat<m,n,S>& A, const mat<n,p,T>& B)
{
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


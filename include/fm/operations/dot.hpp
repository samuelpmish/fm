#pragma once

namespace fm {

template < Kind kind, u32 m, u32 n, typename S, typename T>
constexpr auto dot(const vec<m, S> & x, const matrix<kind, m, n, T> & A) {

  if constexpr (kind == Kind::Isotropic || kind == Kind::Diagonal) {
    static_assert(m == n);
    return x * A.data;
  }

  if constexpr (kind == Kind::Skew) {
    return cross(x, A.data); // ?
  }

  if constexpr (kind == Kind::Symmetric || 
                kind == Kind::Rotation ||
                kind == Kind::Skew ||
                kind == Kind::General) {
    vec<n, decltype(S{} + T{})> output{};
    for (u32 j = 0; j < n; j++) {
      for (u32 i = 0; i < m; i++) {
        output(j) += x(i) * A(i,j);
      }
    }
    return output;
  }

}

template < Kind kind, u32 m, u32 n, typename S, typename T>
constexpr auto dot(const matrix<kind, m, n, S> & A, const vec<n, T> & x) {

  if constexpr (kind == Kind::Isotropic || kind == Kind::Diagonal) {
    static_assert(m == n);
    return A.data * x;
  }

  if constexpr (kind == Kind::Skew) {
    return cross(A.data, x); // ?
  }

  if constexpr (kind == Kind::Symmetric || 
                kind == Kind::Rotation ||
                kind == Kind::Skew ||
                kind == Kind::General) {
    vec<m, decltype(S{} + T{})> output{};
    for (u32 i = 0; i < m; i++) {
      for (u32 j = 0; j < n; j++) {
        output(i) += A(i,j) * x(j);
      }
    }
    return output;
  }

}

template < Kind kindA, Kind kindB, u32 m, u32 n, u32 p, typename TA, typename TB>
constexpr auto dot(const matrix<kindA, m, n, TA> & A, 
                   const matrix<kindB, n, p, TB> & B) {

  using T = decltype(TA{} - TB{});

  if constexpr (kindA == Kind::Isotropic) {
    static_assert(m == p);
    return A.data * B;
  }

  if constexpr (kindA == Kind::Diagonal) {

    if constexpr (kindB == Kind::Isotropic) {
      return A * B.data;
    }

    if constexpr (kindB == Kind::Diagonal) {
      matrix<Kind::Diagonal, m, p, T> output;
      for (u32 i = 0; i < m; i++) {
        output.data[i] = A.data[i] * B.data[i];
      }
      return output;
    }

    if constexpr (kindB == Kind::Symmetric || 
                  kindB == Kind::Rotation ||
                  kindB == Kind::Skew ||
                  kindB == Kind::General) {
      matrix<Kind::General, m, p, T> output;
      for (u32 i = 0; i < m; i++) {
        for (u32 j = 0; j < m; j++) {
          output(i,j) = A.data[i] * B(i,j);
        }
      }
      return output;
    }



  }

}

}
#pragma once

namespace fm {

template < Kind kind, u32 m, u32 n, typename T>
__host__ __device__ constexpr auto dot(const matrix<kind, m, n, T> & A, const double x) {
  return A * x;
}

template < Kind kind, u32 m, u32 n, typename T>
__host__ __device__ constexpr auto dot(const double x, const matrix<kind, m, n, T> & A) {
  return x * A;
}

template < Kind kind, u32 m, u32 n, typename S, typename T>
__host__ __device__ constexpr auto dot(const vec<m, S> & x, const matrix<kind, m, n, T> & A) {

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
    vec<n, decltype(S{} * T{})> output{};
    for (u32 j = 0; j < n; j++) {
      for (u32 i = 0; i < m; i++) {
        output(j) += x(i) * A(i,j);
      }
    }
    return output;
  }

}

template < Kind kind, u32 m, u32 n, typename S, typename T>
__host__ __device__ constexpr auto dot(const matrix<kind, m, n, S> & A, const vec<n, T> & x) {

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
    vec<m, decltype(S{} * T{})> output{};
    for (u32 i = 0; i < m; i++) {
      for (u32 j = 0; j < n; j++) {
        output(i) += A(i,j) * x(j);
      }
    }
    return output;
  }

}

template < Kind kindA, Kind kindB, u32 m, u32 n, u32 p, typename TA, typename TB>
__host__ __device__ constexpr auto dot(const matrix<kindA, m, n, TA> & A, 
                   const matrix<kindB, n, p, TB> & B) {

  using T = decltype(TA{} * TB{});

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
        for (u32 j = 0; j < p; j++) {
          output(i,j) = A.data[i] * B(i,j);
        }
      }
      return output;
    }

  }

  // TODO: specialize remaining cases individually
  if constexpr (kindA == Kind::Symmetric || 
                kindA == Kind::Rotation ||
                kindA == Kind::Skew ||
                kindA == Kind::General) {
    matrix<Kind::General, m, p, T> output;
    for (u32 i = 0; i < m; i++) {
      for (u32 j = 0; j < p; j++) {
        T sum{};
        for (u32 k = 0; k < n; k++) {
          sum += A(i,k) * B(k,j);
        }
        output(i,j) = sum;
      }
    }
    return output;
  }

}

/// returns A(i,j) B(i,j)
template < Kind kindA, Kind kindB, u32 m, u32 n, typename TA, typename TB>
__host__ __device__ constexpr auto ddot(const matrix<kindA, m, n, TA> & A, 
                                        const matrix<kindB, m, n, TB> & B) {

  using T = decltype(TA{} * TB{});

  if constexpr (kindA == Kind::Isotropic) {
    if constexpr (kindB == Kind::Isotropic) {
      return m * A.data * B.data;
    }

    if constexpr (kindB == Kind::Diagonal) {
      T total{};
      for (u32 i = 0; i < m; i++) {
        total += A.data * B.data[i];
      }
      return total;
    }

    if constexpr (kindB == Kind::Skew) {
      return T{}; // identically zero
    }

    if constexpr (kindB == Kind::Rotation) {
      if constexpr (m == 2) { return 2 * A.data * B.c; }
      if constexpr (m == 3) { return A.data * (3.0 - 2.0 * squared_norm(B.s)); }
    }

    if constexpr (kindB == Kind::General || kindB == Kind::Symmetric) {
      T total{};
      for (u32 i = 0; i < m; i++) {
        total += A.data * B(i,i);
      }
      return total;
    }

  }

  if constexpr (kindA == Kind::Diagonal) {
    if constexpr (kindB == Kind::Isotropic) {
      return m * A.data * B.data;
    }

    if constexpr (kindB == Kind::Diagonal) {
      return dot(A.data, B.data);
    }

    if constexpr (kindB == Kind::Skew) {
      return T{}; // identically zero
    }

    if constexpr (kindB == Kind::General || 
                  kindB == Kind::Symmetric || 
                  kindB == Kind::Rotation) {
      T total{};
      for (u32 i = 0; i < m; i++) {
        total += A.data * B(i,i);
      }
      return total;
    }
  }

  if constexpr (kindA == Kind::Skew) {
    if constexpr (kindB == Kind::Isotropic ||
                  kindB == Kind::Diagonal ||
                  kindB == Kind::Symmetric ) {
      return T{}; // identically zero
    }

    if constexpr (kindB == Kind::Skew) {
      return dot(A.data, B.data);
    }

    if constexpr (kindB == Kind::Rotation) {
      if constexpr (m == 2) { return 2.0 * A.data * B.s; }
      if constexpr (m == 3) { return 4.0 * B.c * dot(A.data, B.s); }
    }

    if constexpr (kindB == Kind::General) {
      if constexpr (m == 2) { 
        return A.data * (B(0,1) - B(1,0)); 
      }
      if constexpr (m == 3) { 
        return A.data[0] * (B(2,1) - B(1,2)) + 
               A.data[1] * (B(0,2) - B(2,0)) + 
               A.data[2] * (B(1,0) - B(0,1));
      } 
    }
  }

  // TODO: specialize remaining cases individually
  if constexpr (kindA == Kind::Rotation || 
                kindA == Kind::Symmetric || 
                kindA == Kind::General) {
      T total{};
      for (u32 i = 0; i < m; i++) {
        for (u32 j = 0; j < n; j++) {
          total += A(i,j) * B(i,j);
        }
      }
      return total;
  }

}

/// returns A(i,j) B(i,j)
template < Kind kindA, Kind kindB, u32 m, u32 n, typename TA, typename TB>
__host__ __device__ constexpr auto inner(const matrix<kindA, m, n, TA> & A, 
                                         const matrix<kindB, m, n, TB> & B) {
  return ddot(A, B);
}

}
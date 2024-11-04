#pragma once

namespace fm {

////////////////////////////////////////////////////////////////////////////////

template < Kind kind, u32 rows, u32 cols, typename T>
__host__ __device__ constexpr auto operator*(T a, const matrix<kind, rows, cols, T> & B) {

  if constexpr (kind == Kind::General || kind == Kind::Rotation) {
    mat<rows, cols, T> output;
    for (u32 i = 0; i < rows; i++) {
      for (u32 j = 0; j < cols; j++) {
        output(i,j) = a * B(i,j); 
      }
    }
    return output;
  }

  if constexpr (kind == Kind::Isotropic || 
                kind == Kind::Diagonal || 
                kind == Kind::Skew || 
                kind == Kind::Symmetric) {
    return matrix<kind, rows, cols, T>{a * B.data};
  }

}

template < Kind kind, u32 rows, u32 cols, typename T>
__host__ __device__ constexpr auto operator*(const matrix<kind, rows, cols, T> & B, const T a) {

  if constexpr (kind == Kind::General || kind == Kind::Rotation) {
    mat<rows, cols, T> output;
    for (u32 i = 0; i < rows; i++) {
      for (u32 j = 0; j < cols; j++) {
        output(i,j) = B(i,j) * a;
      }
    }
    return output;
  }

  if constexpr (kind == Kind::Isotropic || 
                kind == Kind::Diagonal || 
                kind == Kind::Skew || 
                kind == Kind::Symmetric) {
    return matrix<kind, rows, cols, T>{B.data * a};
  }

}

template < Kind kind, u32 rows, u32 cols, typename T>
__host__ __device__ constexpr auto operator*=(matrix<kind, rows, cols, T> & A, T scale) {
  static_assert(kind != Kind::Rotation, "rotation matrices don't support operator*=(scalar)");
  if constexpr (kind == Kind::General) {
    for (u32 i = 0; i < rows; i++) {
      for (u32 j = 0; j < cols; j++) {
        A(i,j) *= scale;
      }
    }
  } else {
    A.data = A.data * scale;
  }
}

////////////////////////////////////////////////////////////////////////////////

template < Kind kindA, Kind kindB, u32 m, u32 n, typename T>
__host__ __device__ constexpr auto operator+=(matrix<kindA, m, n, T> & A, 
                          const matrix<kindB, m, n, T> & B) {

  if constexpr (kindA == Kind::General) {
    if constexpr (kindB == Kind::General || 
                  kindB == Kind::Symmetric || 
                  kindB == Kind::Skew || 
                  kindB == Kind::Rotation) {
      for (uint32_t i = 0; i < m; i++) {
        for (uint32_t j = 0; j < n; j++) {
          A(i,j) += B(i,j);
        }
      }
    }

    if constexpr (kindA == Kind::Diagonal) {
      for (uint32_t i = 0; i < m; i++) {
        A(i,i) += B.data[i];
      }
    }

    if constexpr (kindA == Kind::Isotropic) {
      for (uint32_t i = 0; i < m; i++) {
        A(i,i) += B.data;
      }
    }
  }

  if constexpr (kindA == Kind::Symmetric) {
    if constexpr (kindB == Kind::Symmetric) {
      for (uint32_t i = 0; i < decltype(A)::num_values; i++) {
        A.data[i] += B.data[i];
      }
    }

    if constexpr (kindA == Kind::Diagonal) {
      for (uint32_t i = 0; i < m; i++) {
        A(i,i) += B.data[i];
      }
    }

    if constexpr (kindA == Kind::Isotropic) {
      for (uint32_t i = 0; i < m; i++) {
        A(i,i) += B.data;
      }
    }

    if constexpr (kindB == Kind::General || 
                  kindB == Kind::Skew || 
                  kindB == Kind::Rotation) {
      static_assert(always_false<T>{}, "sym += {skew, rotation, general} not supported");
    }
  }

  if constexpr (kindA == Kind::Rotation) {
    static_assert(always_false<T>{}, "rot does not support operator+=");
  }

  if constexpr (kindA == Kind::Skew) {
    if constexpr (kindB == Kind::Skew) {
      for (uint32_t i = 0; i < decltype(A)::num_values; i++) {
        A.data[i] += B.data[i];
      }
    }

    if constexpr (kindB == Kind::General || 
                  kindB == Kind::Symmetric ||
                  kindB == Kind::Diagonal ||
                  kindB == Kind::Isotropic ||
                  kindB == Kind::Rotation) {
      static_assert(always_false<T>{}, "skew += skew is the only supported operator+=()");
    }
  }

  if constexpr (kindA == Kind::Diagonal) {
    if constexpr (kindB == Kind::Diagonal) {
      for (uint32_t i = 0; i < m; i++) {
        A.data[i] += B.data[i];
      }
    }

    if constexpr (kindB == Kind::Isotropic) {
      for (uint32_t i = 0; i < m; i++) {
        A.data[i] += B.data;
      }
    }

    if constexpr (kindB == Kind::General || 
                  kindB == Kind::Symmetric ||
                  kindB == Kind::Skew ||
                  kindB == Kind::Rotation) {
      static_assert(always_false<T>{}, "diag += {mat, sym, skew, rot} are unsupported");
    }
  }

  if constexpr (kindA == Kind::Isotropic) {
    if constexpr (kindB == Kind::Isotropic) {
      A.data += B.data;
    }

    if constexpr (kindB == Kind::General || 
                  kindB == Kind::Symmetric ||
                  kindB == Kind::Skew ||
                  kindB == Kind::Diagonal ||
                  kindB == Kind::Rotation) {
      static_assert(always_false<T>{}, "iso += {diag, mat, sym, skew, rot} are unsupported");
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

template < Kind kind, u32 rows, u32 cols, typename T>
__host__ __device__ constexpr auto operator/(T a, const matrix<kind, rows, cols, T> & B) {

  if constexpr (kind == Kind::General || kind == Kind::Rotation) {
    mat<rows, cols, T> output;
    for (u32 i = 0; i < rows; i++) {
      for (u32 j = 0; j < cols; j++) {
        output(i,j) = a / B(i,j); 
      }
    }
    return output;
  }

  if constexpr (kind == Kind::Isotropic || 
                kind == Kind::Diagonal || 
                kind == Kind::Skew || 
                kind == Kind::Symmetric) {
    return matrix<kind, rows, cols, T>{a / B.data};
  }

}

template < Kind kind, u32 rows, u32 cols, typename T>
__host__ __device__ constexpr auto operator/(const matrix<kind, rows, cols, T> & B, const T a) {

  if constexpr (kind == Kind::General || kind == Kind::Rotation) {
    mat<rows, cols, T> output;
    for (u32 i = 0; i < rows; i++) {
      for (u32 j = 0; j < cols; j++) {
        output(i,j) = B(i,j) / a;
      }
    }
    return output;
  }

  if constexpr (kind == Kind::Isotropic || 
                kind == Kind::Diagonal || 
                kind == Kind::Skew || 
                kind == Kind::Symmetric) {
    return matrix<kind, rows, cols, T>{B.data / a};
  }

}

template < Kind kind, u32 rows, u32 cols, typename T>
__host__ __device__ constexpr auto operator/=(matrix<kind, rows, cols, T> & A, T scale) {
  static_assert(kind != Kind::Rotation, "rotation matrices don't support operator/=(scalar)");
  if constexpr (kind == Kind::General) {
    for (u32 i = 0; i < rows; i++) {
      for (u32 j = 0; j < cols; j++) {
        A(i,j) /= scale;
      }
    }
  } else {
    A.data = A.data / scale;
  }
}

////////////////////////////////////////////////////////////////////////////////

template < Kind kindA, Kind kindB, u32 rows, u32 cols, typename TA, typename TB>
__host__ __device__ constexpr auto operator-(const matrix<kindA, rows, cols, TA> & A, 
                         const matrix<kindB, rows, cols, TB> & B) {

  using T = decltype(TA{} - TB{});

  if constexpr (kindA == Kind::General || kindB == Kind::General) {
    mat<rows, cols, T> output;
    for (u32 i = 0; i < rows; i++) {
      for (u32 j = 0; j < cols; j++) {
        output(i,j) = A(i,j) - B(i,j); 
      }
    }
    return output;
  }

  if constexpr (kindA == Kind::Symmetric) {
    if constexpr (kindB == Kind::Symmetric) {
      return matrix< Kind::Symmetric, rows, cols, T >{A.data - B.data};
    }

    if constexpr (kindB == Kind::Isotropic) {
      matrix< Kind::Symmetric, rows, cols, T > output{};
      for (u32 i = 0; i < decltype(output)::num_values; i++) {
        output.data[i] = A.data[i];
      }

      for (u32 i = 0; i < rows; i++) {
        output(i, i) -= B.data;
      }
      return output;
    }
  }

  if constexpr (kindA == Kind::Isotropic && kindB == Kind::Isotropic) {
    return iso<rows, T>{A.data - B.data};
  }

  if constexpr (kindA == Kind::Diagonal && kindB == Kind::Diagonal) {
    return diag<rows, T>{A.data - B.data};
  }

  if constexpr (kindA == Kind::Skew && kindB == Kind::Skew) {
    return skew<rows, T>{A.data - B.data};
  }

}

template < Kind kindA, Kind kindB, u32 rows, u32 cols, typename TA, typename TB>
__host__ __device__ constexpr auto operator+(const matrix<kindA, rows, cols, TA> & A, 
                         const matrix<kindB, rows, cols, TB> & B) {

  using T = decltype(TA{} + TB{});

  if constexpr (kindA == Kind::General || kindB == Kind::General) {
    mat<rows, cols, T> output;
    for (u32 i = 0; i < rows; i++) {
      for (u32 j = 0; j < cols; j++) {
        output(i,j) = A(i,j) + B(i,j); 
      }
    }
    return output;
  }

  if constexpr (kindA == Kind::Isotropic && kindB == Kind::Isotropic) {
    return iso<rows, T>{A.data + B.data};
  }

  if constexpr (kindA == Kind::Diagonal && kindB == Kind::Diagonal) {
    return diag<rows, T>{A.data + B.data};
  }

  if constexpr (kindA == Kind::Skew && kindB == Kind::Skew) {
    return skew<rows, T>{A.data + B.data};
  }

}

}
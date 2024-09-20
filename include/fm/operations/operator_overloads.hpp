#pragma once

namespace fm {

////////////////////////////////////////////////////////////////////////////////

template < Kind kind, u32 rows, u32 cols, typename T>
constexpr auto operator*(T a, const matrix<kind, rows, cols, T> & B) {

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
constexpr auto operator*(const matrix<kind, rows, cols, T> & B, const T a) {

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
constexpr auto operator*=(matrix<kind, rows, cols, T> & A, T scale) {
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

template < Kind kind, u32 rows, u32 cols, typename T>
constexpr auto operator/(T a, const matrix<kind, rows, cols, T> & B) {

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
constexpr auto operator/(const matrix<kind, rows, cols, T> & B, const T a) {

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
constexpr auto operator/=(matrix<kind, rows, cols, T> & A, T scale) {
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
constexpr auto operator-(const matrix<kindA, rows, cols, TA> & A, 
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
constexpr auto operator+(const matrix<kindA, rows, cols, TA> & A, 
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
    return iso<rows, T>{A.data - B.data};
  }

  if constexpr (kindA == Kind::Diagonal && kindB == Kind::Diagonal) {
    return diag<rows, T>{A.data - B.data};
  }

  if constexpr (kindA == Kind::Skew && kindB == Kind::Skew) {
    return skew<rows, T>{A.data - B.data};
  }

}

}
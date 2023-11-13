#pragma once

namespace fm {

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

}
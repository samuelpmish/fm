#pragma once

template < uint32_t n, typename T = double >
struct sym { 
  using data_type = T;
  static constexpr int num_values = (n*(n+1)) / 2;
  static constexpr int dimensions[2] = {n,n};
  static constexpr fm::type type = fm::type::sym;

  T data[num_values];

  auto & operator()(uint32_t i, uint32_t j) { return data[index(i,j)]; }
  const auto & operator()(uint32_t i, uint32_t j) const { return data[index(i,j)]; }
  uint32_t index(uint32_t i, uint32_t j) const {
    uint32_t i_upper = std::min(i, j);
    uint32_t j_upper = std::max(i, j);
    return j_upper + ((2 * n - i_upper - 1) * i_upper) / 2;
  }
};

using sym2 = sym<2>;
using sym3 = sym<3>;

template < uint32_t dim, typename T >
constexpr auto operator*(const double & u, const sym< dim, T > & v) {
  sym< dim, T > out;
  for (uint32_t i = 0; i < out.num_values; i++) {
    out.data[i] = u * v.data[i];
  }
  return out; 
}

template < uint32_t dim, typename T >
constexpr auto operator*(const sym< dim, T > & v, const double & u) {
  sym< dim, T > out{};
  for (uint32_t i = 0; i < out.num_values; i++) {
    out.data[i] = u * v.data[i];
  }
  return out; 
}

//template < uint32_t dim, typename T >
//constexpr auto dev(const sym< dim, T > & A) {
//  return A - (tr(A) / dim) * I;
//}
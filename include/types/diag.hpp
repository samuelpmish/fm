#pragma once

template < uint32_t dim, typename T = float >
struct diag { 
  using data_type = T;
  static constexpr int dimensions[2] = {dim,dim};
  static constexpr fm::type type = fm::type::diag;

  vec< dim, T > data; 
};

template < uint32_t dim, typename T = float >
diag(vec<dim, T>) -> diag<dim, T>;

template < uint32_t dim, typename S, typename T >
auto dot(const diag<dim, S> & A, const diag<dim, T> & B) {
  return diag{A.data * B.data};
}

template < uint32_t dim, typename S, typename T >
auto operator+(const diag<dim, S> & A, const diag<dim, T> & B) {
  return diag{A.data + B.data};
}

template < uint32_t dim, typename S, typename T >
auto operator-(const diag<dim, S> & A, const diag<dim, T> & B) {
  return diag{A.data - B.data};
}

template < uint32_t dim, typename S, typename T >
auto operator*(const diag<dim, S> & A, const diag<dim, T> & B) {
  return diag{A.data * B.data};
}

template < uint32_t dim, typename S, typename T >
auto operator/(const diag<dim, S> & A, const diag<dim, T> & B) {
  return diag{A.data / B.data};
}

template < uint32_t dim, typename T >
auto inv(const diag<dim, T> & A) {
  return diag<dim>{T(1) / A.data};
}

template < uint32_t dim, typename T >
auto det(const diag<dim, T> & A) {
  T out = 1;
  for (uint32_t i = 0; i < dim; i++) { out *= A.data[i]; }
  return out;
}
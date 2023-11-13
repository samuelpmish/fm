#pragma once

#include <cmath>

template < uint32_t dim, typename T = double >
struct iso { 
  using data_type = T;
  static constexpr int dimensions[2] = {dim,dim};
  static constexpr fm::type type = fm::type::iso;

  T data; 
};

template < uint32_t dim, typename S, typename T >
auto dot(const iso<dim, S> & A, const iso<dim, T> & B) {
  return iso<dim>{A.data * B.data};
}

template < uint32_t dim, typename S, typename T >
auto operator+(const iso<dim, S> & A, const iso<dim, T> & B) {
  return iso<dim>{A.data + B.data};
}

template < uint32_t dim, typename S, typename T >
auto operator-(const iso<dim, S> & A, const iso<dim, T> & B) {
  return iso<dim>{A.data - B.data};
}

template < uint32_t dim, typename S, typename T >
auto operator*(const iso<dim, S> & A, const iso<dim, T> & B) {
  return iso<dim>{A.data * B.data};
}

template < uint32_t dim, typename T >
auto operator*(double scale, const iso<dim, T> & A) {
  return iso<dim, T>{scale * A.data};
}

template < uint32_t dim, typename S, typename T >
auto operator/(const iso<dim, S> & A, const iso<dim, T> & B) {
  return iso<dim>{A.data / B.data};
}

template < uint32_t dim, typename T >
auto inv(const iso<dim, T> & A) {
  return iso<dim>{T(1) / A.data};
}

template < uint32_t dim, typename T >
auto det(const iso<dim, T> & A) {
  return std::pow(A.data, dim);
}
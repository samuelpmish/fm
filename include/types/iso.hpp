#pragma once

#include <cmath>

template < uint32_t dim, typename T = float >
struct iso { T data; };

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
#pragma once

template < uint32_t dim, typename T = double >
struct ortho;

template < typename T >
struct ortho< 2, T > {
  using data_type = T;
  static constexpr int dimensions[2] = {2, 2};
  static constexpr fm::type type = fm::type::ortho;

  T c;
  T s;
};

template < typename T >
struct ortho< 3, T > {
  using data_type = T;
  static constexpr int dimensions[2] = {3, 3};
  static constexpr fm::type type = fm::type::ortho;

  T c;
  vec< 3, T > s;
};

using ortho2 = ortho<2, double>;
using ortho3 = ortho<3, double>;

using ortho2f = ortho<2, float>;
using ortho3f = ortho<3, float>;

template < typename T >
ortho< 2, T > rotation_matrix(T theta) {
  return ortho< 2, T >{std::cos(theta), std::sin(theta)};
};

template < typename T >
ortho< 3, T > rotation_matrix(vec<3,T> axis) {
  T theta = norm(axis);
  return ortho< 3, T >{std::cos(theta / 2), std::sin(theta / 2) * normalize(axis)};
};

template < typename T >
mat< 2, 2, T > as_mat(const ortho<2, T> & R) {
  auto [c, s] = R;
  return mat< 2, 2, T >{{{c, -s}, {s, c}}};
}

template < typename T >
mat< 3, 3, T > as_mat(const ortho<3, T> & R) {
  auto [c, s] = R;
  return mat< 3, 3, T >{{
    {1-2*(s[1]*s[1]+s[2]*s[2]),   2*(s[0]*s[1]-   c*s[2]),   2*(s[0]*s[2]+   c*s[1])},
    {  2*(s[0]*s[1]+   c*s[2]), 1-2*(s[0]*s[0]+s[2]*s[2]),   2*(s[1]*s[2]-   c*s[0])},
    {  2*(s[0]*s[2]-   c*s[1]),   2*(s[1]*s[2]+   c*s[0]), 1-2*(s[0]*s[0]+s[1]*s[1])}
  }};
}

template < uint32_t dim, typename T >
ortho<dim,T> transpose(const ortho<dim,T> & Q) { return ortho<dim,T>{Q.c, -Q.s}; }

template < uint32_t dim, typename T >
ortho<dim,T> inv(const ortho<dim,T> & Q) { return transpose(Q); }
#pragma once

#include "fm/types/vec.hpp"

#include <algorithm>

namespace fm {

template < int dim, typename T = float >
struct AABB { 
  vec<dim,T> min; 
  vec<dim,T> max; 

  __host__ __device__ T SDF(const vec<dim,T> & p) const {
    constexpr T zero{};
    vec<dim,T> center = (min + max) * 0.5f;
    vec<dim,T> halfwidths = (max - min) * 0.5f;
    vec<dim,T> q = abs(p - center) - halfwidths;
    return norm(fm::max(q, zero)) + std::min(fm::max(q), zero);
  }
};

template < int dim, typename T >
__host__ __device__ AABB<dim,T> bounding_box(AABB<dim,T> a, AABB<dim,T> b) {
  return AABB<dim,T>{fm::min(a.min, b.min), fm::max(a.max, b.max)};
}

template < typename T >
__host__ __device__ bool intersecting(AABB<2,T> a, AABB<2,T> b) {
  return a.min[0] <= b.max[0] && a.min[1] <= b.max[1] && a.max[0] >= b.min[0] && a.max[1] >= b.min[1];
}

template < typename T >
__host__ __device__ bool intersecting(AABB<3,T> a, AABB<3,T> b) {
  return a.min[0] <= b.max[0] && a.min[1] <= b.max[1] && a.min[2] <= b.max[2] && 
         a.max[0] >= b.min[0] && a.max[1] >= b.min[1] && a.max[2] >= b.min[2];
}

template < int dim, typename T >
__host__ __device__ AABB<dim,T> intersection_of(AABB<dim,T> a, AABB<dim,T> b) {
  return AABB<dim,T>{fm::max(a.min, b.min), fm::min(a.max, b.max)};
}

}
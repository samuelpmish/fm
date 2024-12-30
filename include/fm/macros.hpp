#pragma once

#ifndef __CUDACC__
#define __host__
#define __device__
#define __constant__
#endif

namespace fm {
  template < typename T >
  __host__ __device__ void swap(T & A, T & B) {
    T tmp = B;
    B = A;
    A = tmp;
  }
}
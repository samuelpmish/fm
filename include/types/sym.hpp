#pragma once

template < uint32_t n, typename T = float >
struct sym { 
  T data[(n*(n+1)) / 2];

  auto & operator[](uint32_t i) { return data[i]; }
  const auto & operator[](uint32_t i) const { return data[i]; }
};
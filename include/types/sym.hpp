#pragma once

template < uint32_t n, typename T = float >
struct sym { 
  T data[(n*(n+1)) / 2];

  auto & operator()(uint32_t i, uint32_t j) { return data[index(i,j)]; }
  const auto & operator()(uint32_t i, uint32_t j) const { return data[index(i,j)]; }
  uint32_t index(uint32_t i, uint32_t j) {
    uint32_t i_upper = std::min(i, j);
    uint32_t j_upper = std::max(i, j);
    return j_upper + ((2 * n - i_upper - 1) * i_upper) / 2;
  }
};

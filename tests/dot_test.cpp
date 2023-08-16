#include "common.hpp"

#include "operations/dot.hpp"

TEST(UnitTest, dot) {

  vec2 u = {1.0, 2.0};

  sym<2, double> A{{
    1.0, 2.0,
         3.0
  }};

  mat<2,3,double> B{{
    {1.0, 2.0, 3.0},
    {4.0, 5.0, 6.0}
  }};

  vec3 v = {1.0, 2.0, 3.0};

  vec2 uA = dot(u, A); 
  EXPECT_EQ(uA, (vec2{5.0, 8.0}));

  vec3 uAB = dot(uA, B); 
  EXPECT_EQ(uAB, (vec3{37.0, 50.0, 63.0}));

  double uABv = dot(uAB, v);
  EXPECT_EQ(uABv, 326.0);

  vec2 Bv = dot(B, v);
  EXPECT_EQ(Bv, (vec2{14.0, 32.0}));

  vec2 ABv = dot(A, Bv);
  EXPECT_EQ(ABv, (vec2{78.0, 124.0}));

  uABv = dot(u, ABv);
  EXPECT_EQ(uABv, 326.0);

  mat2 AA = dot(A, A);
  EXPECT_EQ(AA, (mat2{{{5.0, 8.0}, {8.0, 13.0}}}));

  mat2 BBT = dot(B, transpose(B));
  EXPECT_EQ(BBT, (mat2{{{14., 32.}, {32., 77.}}}));

  mat3 BTB = dot(transpose(B), B);
  EXPECT_EQ(BTB, (mat3{{{17., 22., 27.}, {22., 29., 36.}, {27., 36., 45.}}}));

}

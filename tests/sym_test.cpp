#include "common.hpp"

TEST(UnitTest, sym2) {
  sym<2, double> A2{{
    1.0, 2.0,
         3.0
  }};
  EXPECT_EQ(A2.index(0, 0), 0);
  EXPECT_EQ(A2.index(0, 1), 1);
  EXPECT_EQ(A2.index(1, 1), 2);

  EXPECT_EQ(A2.data[0], 1);
  EXPECT_EQ(A2.data[1], 2);
  EXPECT_EQ(A2.data[2], 3);
}

TEST(UnitTest, sym3) {
  sym<3, double> A3{{
    1.0, 2.0, 3.0,
         4.0, 5.0,
              6.0
  }};

  EXPECT_EQ(A3.index(0, 0), 0);
  EXPECT_EQ(A3.index(0, 1), 1);
  EXPECT_EQ(A3.index(0, 2), 2);
  EXPECT_EQ(A3.index(1, 1), 3);
  EXPECT_EQ(A3.index(1, 2), 4);
  EXPECT_EQ(A3.index(2, 1), 4);
  EXPECT_EQ(A3.index(2, 2), 5);

  for (uint32_t i = 0; i < 6; i++) {
    EXPECT_EQ(A3.data[i], i+1);
  }

  EXPECT_EQ(A3(0,0) + A3(1,0) + A3(2, 0), 6);

}

TEST(UnitTest, sym6) {
  sym<6, double> A6 = {{
    1.0, 2.0,  3.0,  4.0,  5.0,  6.0,
         7.0,  8.0,  9.0, 10.0, 11.0,
              12.0, 13.0, 14.0, 15.0,
                    16.0, 17.0, 18.0,
                          19.0, 20.0,
                                21.0
  }};
  EXPECT_EQ(A6.index(0, 0),  0);
  EXPECT_EQ(A6.index(0, 1),  1);
  EXPECT_EQ(A6.index(0, 2),  2);
  EXPECT_EQ(A6.index(0, 3),  3);
  EXPECT_EQ(A6.index(0, 4),  4);
  EXPECT_EQ(A6.index(0, 5),  5);

  EXPECT_EQ(A6.index(1, 1),  6);
  EXPECT_EQ(A6.index(1, 2),  7);
  EXPECT_EQ(A6.index(1, 3),  8);
  EXPECT_EQ(A6.index(1, 4),  9);
  EXPECT_EQ(A6.index(1, 5), 10);

  EXPECT_EQ(A6.index(2, 2), 11);
  EXPECT_EQ(A6.index(2, 3), 12);
  EXPECT_EQ(A6.index(2, 4), 13);
  EXPECT_EQ(A6.index(2, 5), 14);

  EXPECT_EQ(A6.index(3, 3), 15);
  EXPECT_EQ(A6.index(3, 4), 16);
  EXPECT_EQ(A6.index(3, 5), 17);

  EXPECT_EQ(A6.index(4, 4), 18);
  EXPECT_EQ(A6.index(4, 5), 19);

  EXPECT_EQ(A6.index(5, 5), 20);

  for (uint32_t i = 0; i < 21; i++) {
    EXPECT_EQ(A6.data[i], i+1);
  }

  EXPECT_EQ(A6(0,0) + A6(1,0) + A6(2, 0), 6);

}

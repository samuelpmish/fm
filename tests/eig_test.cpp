#include "common.hpp"

#include "fm/operations/eig.hpp"

TEST(eig, identity) {
  sym3 A{{
    1.0, 0.0, 0.0,
         1.0, 0.0,
              1.0
  }};

  auto [evalues, evectors] = eig(A);

  mat3 B = dot(dot(transpose(evectors), diag3{evalues}), evectors);

  mat3 error = A - B;

  EXPECT_NEAR(norm(error), 0.0, 1.0e-14);
}

TEST(eig, laplace) {
  sym3 A{{
   -2.0,  1.0,  0.0,
         -2.0,  1.0,
               -2.0
  }};

  auto [evalues, evectors] = eig(A);

  mat3 B = dot(dot(transpose(evectors), diag3{evalues}), evectors);

  mat3 error = A - B;

  EXPECT_NEAR(norm(error), 0.0, 1.0e-14);
}

TEST(eig, repeated_eigenvalues) {
  sym3 A{{
    2.0,  0.0,  0.0,
          2.0,  0.0,
                1.0
  }};

  auto [evalues, evectors] = eig(A);

  mat3 B = dot(dot(transpose(evectors), diag3{evalues}), evectors);

  mat3 error = A - B;

  EXPECT_NEAR(norm(error), 0.0, 1.0e-14);
}

TEST(eig, nearly_repeated_eigenvalues) {
  sym3 A{{
    2.0,  1.0e-9,  0.0,
          2.0,      0.0,
                    1.0
  }};

  std::cout << "nearly repeated" << std::endl;

  auto [evalues, evectors] = eig(A);

  mat3 B = dot(dot(transpose(evectors), diag3{evalues}), evectors);

  std::cout << std::setprecision(16) << evalues << std::endl;

  mat3 error = A - B;

  EXPECT_NEAR(norm(error) / norm(A), 0.0, 1.0e-12);
}

TEST(eig, rank1_matrix) {
  sym3 A{{
    1.0, 1.0, 1.0,
         1.0, 1.0,
              1.0
  }};

  auto [evalues, evectors] = eig(A);

  mat3 B = dot(dot(transpose(evectors), diag3{evalues}), evectors);

  mat3 error = A - B;

  EXPECT_NEAR(norm(error) / norm(A), 0.0, 1.0e-13);
}

TEST(eig, rank2_matrix) {
  sym3 A{{
    2.0, 1.0, 0.0,
         2.0, 0.0,
              0.0
  }};

  auto [evalues, evectors] = eig(A);

  mat3 B = dot(dot(transpose(evectors), diag3{evalues}), evectors);

  mat3 error = A - B;

  EXPECT_NEAR(norm(error) / norm(A), 0.0, 1.0e-13);
}

TEST(eig, offdiagonal) {
  sym3 A{{
    0.0, 1.0, 2.0,
         0.0, 3.0,
              0.0
  }};

  auto [evalues, evectors] = eig(A);

  mat3 B = dot(dot(transpose(evectors), diag3{evalues}), evectors);

  mat3 error = A - B;

  EXPECT_NEAR(norm(error) / norm(A), 0.0, 1.0e-13);
}

TEST(eig, randomized_tests) {
  
  double eps = 1.0e-9;

  vec3 lambda = {1.0, 2.0 + eps, 2.0 - eps};
  //mat3 R = RotationMatrix(vec3{0.5, 1.7, 2.9}); 
  mat3 R = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};

  std::cout << dot(dot(transpose(R), diag3{lambda}), R) << std::endl;
  sym3 A = sym3(dot(dot(transpose(R), diag3{lambda}), R));
  std::cout << A << std::endl;

  auto [evalues, evectors] = eig(A);

  mat3 B = dot(dot(transpose(evectors), diag3{evalues}), evectors);

  mat3 error = A - B;

  std::cout << A << std::endl;

  std::cout << evalues << std::endl;
  std::cout << evectors << std::endl;
  std::cout << norm(error) / norm(A) << std::endl;

  EXPECT_NEAR(norm(error) / norm(A), 0.0, 1.0e-13);
}
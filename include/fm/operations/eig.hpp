#pragma once

#include "fm/operations/invariants.hpp"
#include "fm/operations/operator_overloads.hpp"

namespace fm {

__host__ __device__ constexpr double signum(double x) {
  if (x < 0) return -1.0;
  if (x > 0) return +1.0;
  return 0.0; // if x == 0.0
}

template < uint32_t n, typename T = double >
struct Eigensystem {
  vec<n, T> lambda;
  mat<n, n, T> X;
};

// eigendecomposition for symmetric A
//
// based on "A robust algorithm for finding the eigenvalues and
// eigenvectors of 3x3 symmetric matrices", by Scherzinger & Dohrmann
__host__ __device__ inline Eigensystem<3> eig(const sym3 & A) {

  vec<3> eigenvalues = {0.0, 0.0, 0.0};
  mat<3,3> eigenvectors = {{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};

  double trA_3 = invariant<1>(A) / 3.0;
  sym3 A_dev = A - trA_3 * iso3{1.0};
  double J2 = invariant<2>(A_dev);
  double J3 = invariant<3>(A_dev);

  if (J2 > 0.0) {

    // angle used to find eigenvalues
    double tmp = (0.5 * J3) * pow(3.0 / J2, 1.5);
    double alpha = acos(fmin(fmax(tmp, -1.0), 1.0)) / 3.0;

    // consider the most distinct eigenvalue first
    if (6.0 * alpha < M_PI) {
      eigenvalues[0] = 2 * sqrt(J2 / 3.0) * cos(alpha);
    } else {
      eigenvalues[0] = 2 * sqrt(J2 / 3.0) * cos(alpha + 2.0 * M_PI / 3.0);
    }

    std::cout << J2 << " " << J3 << std::endl;
    std::cout << eigenvalues[0] << std::endl;

    // find the eigenvector for that eigenvalue
    vec3 r[3];

    int imax = -1;
    double norm_max = -1.0;

    for (int i = 0; i < 3; i++) {

      for (int j = 0; j < 3; j++) {
        r[i][j] = A_dev(j, i) - (i == j) * eigenvalues[0];
      }

      double norm_r = norm(r[i]);
      if (norm_max < norm_r) {
        imax = i;
        norm_max = norm_r;
      }

    }

    vec3 s0, s1, t1, t2, v0, v1, v2, w;

    s0 = normalize(r[imax]);
    t1 = r[(imax+1)%3] - s0 * dot(r[(imax+1)%3], s0);
    t2 = r[(imax+2)%3] - s0 * dot(r[(imax+2)%3], s0);
    s1 = normalize((norm(t1) > norm(t2)) ? t1 : t2);

    // record the first eigenvector
    v0 = cross(s0, s1);
    eigenvectors[0] = v0;

    // get the other two eigenvalues by solving the
    // remaining quadratic characteristic polynomial
    auto A_dev_s0 = dot(A_dev, s0);
    auto A_dev_s1 = dot(A_dev, s1);

    double A11 = dot(s0, A_dev_s0);
    double A12 = dot(s0, A_dev_s1);
    double A21 = dot(s1, A_dev_s0);
    double A22 = dot(s1, A_dev_s1);

#if 1
    //double delta = 0.5 * signum(A11-A22) * sqrt((A11-A22)*(A11-A22) + 4*A12*A21);
    double delta = 0.5 * sqrt((A11-A22)*(A11-A22) + 4*A12*A21);

    std::cout << A11 << " " << A12 << " " << A21 << " " << A22 << " " << delta << std::endl;

    eigenvalues[1] = 0.5 * (A11 + A22) - delta;
    eigenvalues[2] = 0.5 * (A11 + A22) + delta;
#else
    eigenvalues[1] = 0.5 * (A11 + A22) - 0.5 * signum(A11-A22) * sqrt((A11-A22)*(A11-A22) + 4*A12*A21);
    eigenvalues[2] = 0.5 * (A11 + A22) - eigenvalues[1];
#endif

    // if the remaining eigenvalues are exactly the same
    // then just use the basis for the orthogonal complement
    // found earlier
    if (fabs(delta) <= 1.0e-15) {
      
      eigenvectors[1] = s0;
      eigenvectors[2] = s1;
      
    // otherwise compute the remaining eigenvectors
    } else {

      t1 = A_dev_s0 - s0 * eigenvalues[1];
      t2 = A_dev_s1 - s1 * eigenvalues[1];

      w = normalize((norm(t1) > norm(t2)) ? t1 : t2);

      v1 = normalize(cross(w, v0));
      eigenvectors[1] = v1;

      // define the last eigenvector as
      // the direction perpendicular to the
      // first two directions
      v2 = normalize(cross(v0, v1));
      eigenvectors[2] = v2;

    }

  }

  // eta are actually eigenvalues of A_dev, so
  // shift them to get eigenvalues of A
  eigenvalues = eigenvalues + vec3{1,1,1} * trA_3;

  return Eigensystem<3>{ eigenvalues, eigenvectors };
  
}

}
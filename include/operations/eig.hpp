#pragma once

#include "types/sym.hpp"
#include "types/diag.hpp"
#include "types/ortho.hpp"

// eigendecomposition for symmetric A
//
// based on "A robust algorithm for finding the eigenvalues and
// eigenvectors of 3x3 symmetric matrices", by Scherzinger & Dohrmann
inline Eigensystem eig(const sym3 & A) {

  diag<3> eigenvalues = {1.0, 1.0, 1.0};
  mat<3,3> eigenvectors = mat3::Identity();

  auto A_dev = dev(A);

  double J2 = I2(A_dev);
  double J3 = I3(A_dev);

  if (J2 > 0.0) {

    // angle used to find eigenvalues
    double tmp = (0.5 * J3) * pow(3.0 / J2, 1.5);
    double alpha = acos(fmin(fmax(tmp, -1.0), 1.0)) / 3.0;

    // consider the most distinct eigenvalue first
    if (6.0 * alpha < M_PI) {
      eta(0) = 2 * sqrt(J2 / 3.0) * cos(alpha);
    } else {
      eta(0) = 2 * sqrt(J2 / 3.0) * cos(alpha + 2.0 * M_PI / 3.0);
    }

    // find the eigenvector for that eigenvalue
    r1tensor < 3 > r[3];

    int imax = -1;
    double norm_max = -1.0;

    for (int i = 0; i < 3; i++) {

      for (int j = 0; j < 3; j++) {
        r[i](j) = A_dev(j,i) - (i == j) * eta(0);
      }

      double norm_r = norm(r[i]);
      if (norm_max < norm_r) {
        imax = i;
        norm_max = norm_r;
      }

    }

    r1tensor < 3 > s0, s1, t1, t2, v0, v1, v2, w;

    s0 = normalize(r[imax]);
    t1 = r[(imax+1)%3] - dot(r[(imax+1)%3], s0) * s0;
    t2 = r[(imax+2)%3] - dot(r[(imax+2)%3], s0) * s0;
    s1 = normalize((norm(t1) > norm(t2)) ? t1 : t2);

    // record the first eigenvector
    v0 = cross(s0, s1);
    for (int i = 0; i < 3; i++) {
      Q(i,0) = v0(i);
    }

    // get the other two eigenvalues by solving the
    // remaining quadratic characteristic polynomial
    auto A_dev_s0 = dot(A_dev, s0);
    auto A_dev_s1 = dot(A_dev, s1);

    double A11 = dot(s0, A_dev_s0);
    double A12 = dot(s0, A_dev_s1);
    double A21 = dot(s1, A_dev_s0);
    double A22 = dot(s1, A_dev_s1);

    double delta = 0.5 * signum(A11-A22) * sqrt((A11-A22)*(A11-A22) + 4*A12*A21);

    eta(1) = 0.5 * (A11 + A22) - delta;
    eta(2) = 0.5 * (A11 + A22) + delta;

    // if the remaining eigenvalues are exactly the same
    // then just use the basis for the orthogonal complement
    // found earlier
    if (fabs(delta) <= 1.0e-15) {
      
      for (int i = 0; i < 3; i++){
        Q(i,1) = s0(i);
        Q(i,2) = s1(i);
      } 
      
    // otherwise compute the remaining eigenvectors
    } else {

      t1 = A_dev_s0 - eta(1) * s0;
      t2 = A_dev_s1 - eta(1) * s1;

      w = normalize((norm(t1) > norm(t2)) ? t1 : t2);

      v1 = normalize(cross(w, v0));
      for (int i = 0; i < 3; i++) Q(i,1) = v1(i);

      // define the last eigenvector as
      // the direction perpendicular to the
      // first two directions
      v2 = normalize(cross(v0, v1));
      for (int i = 0; i < 3; i++) Q(i,2) = v2(i);

    }

    // eta are actually eigenvalues of A_dev, so
    // shift them to get eigenvalues of A
    eta += tr(A) / 3.0;

  }

}
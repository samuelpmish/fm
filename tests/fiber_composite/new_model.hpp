#pragma once

#include "fm.hpp"
#include "operations/dot.hpp"
#include "operations/linear_solve.hpp"

using vec6 = vec<6, double>;

vec6 voigt(const mat3 & A) {
  return vec6{A(0,0), A(1,1), A(2,2), A(1,2), A(0,2), A(0,1)};
}

mat3 voigt(const vec6 & v) {
  return mat3{{
    {v[0], v[5], v[4]},
    {v[5], v[1], v[3]},
    {v[4], v[3], v[2]}
  }};
}

struct RotatedFiberCompositeModel {

  struct Parameters {
    double Em;
    double Ez;
    double Ep;
    double num;
    double nuzp;
    double nupp;
    double Gzp;
    double Vf;
    double am;
    double afzz;
    double afpp;
    bool thermalStress;
    bool cylDom;

    vec3 compute_thermal_expansion_coefficients() {

      if (thermalStress) {

        double Vm = 1.0 - Vf;

        double s0 = Vf/Ep - (nuzp*nuzp*Vf)/Ez + (Vm - num*num*Vm)/Em;
        double s1 = -((nupp*Vf)/Ep) - (nuzp*nuzp*Vf)/Ez - (num*(1 + num)*Vm)/Em;
        mat2 Stildbar = {{{s0, s1}, {s1, s0}}};

        vec2 Dtildn = {nuzp * Vf + num * Vm, nuzp * Vf + num * Vm};
        vec2 Btildbar = {nuzp * Vf + num * Vm, nuzp * Vf + num * Vm}; // the same as Dtildn?
        vec2 alphatildbar = {afpp * Vf + afzz * nuzp * Vf + am * (1 + num) * Vm, afpp * Vf + afzz * nuzp * Vf + am * (1 + num) * Vm};

        mat2 C22 = inv(Stildbar);
        vec2 C12 = dot(Dtildn, C22);
        vec2 C21 = dot(C22, Btildbar);
        double C11 = Em*Vm + Ez*Vf + dot(Dtildn, C12);

        mat3 Ceff = {{
          {C11,    C12[0],    C12[1]},
          {C21[0], C22[0][0], C22[0][1]},
          {C21[1], C22[1][0], C22[1][1]}
        }};

        double alpha1 = Em*Vm*am + Ez*Vf*afzz + dot(Dtildn, alphatildbar);
        vec2 alpha2 = dot(C22, alphatildbar);

        vec3 alphabar = {alpha1, alpha2[0], alpha2[1]};

        return linear_solve(Ceff, alphabar);

      } else {

        return vec3{};

      }

    }

    mat<6,6,double> compute_stiffness_tensor() {
      auto Vm = 1-Vf; // volume fraction

      double j1 = 1.0/Em + (Vf/Vm)*(1.0/Ep - nuzp*nuzp/Ez);
      double j2 = -num/Em - (Vf/Vm)*(nupp/Ep + nuzp*nuzp/Ez);
      double j3 = 2.0 * ((1.0+num)/Em + (Vf/Vm)*(1.0/(2.0*Gzp)));
      double j4 = 2.0 * ((1.0+num)/Em + (Vf/Vm)*(1.0+nupp)/Ep);

      mat<6,6,double> Jmat = {{
        {       j1,         j2,  -num / Em,  0.0, 0.0, 0.0},
        {       j2,         j1,  -num / Em,  0.0, 0.0, 0.0},
        {-num / Em, - num / Em,   1.0 / Em,  0.0, 0.0, 0.0},
        {      0.0,        0.0,        0.0,   j3, 0.0, 0.0},
        {      0.0,        0.0,        0.0,  0.0,  j3, 0.0},
        {      0.0,        0.0,        0.0,  0.0, 0.0,  j4}
      }};

      mat<6,6,double> Bmat = {{
        {1.0 / Vm,    0.0  , (Vf/Vm)*nuzp,    0.0  ,    0.0  ,    0.0  },
        {     0.0, 1.0 / Vm, (Vf/Vm)*nuzp,    0.0  ,    0.0  ,    0.0  },
        {     0.0,    0.0  ,      1.0    ,    0.0  ,    0.0  ,    0.0  },  
        {     0.0,    0.0  ,      0.0    , 1.0 / Vm,    0.0  ,    0.0  },
        {     0.0,    0.0  ,      0.0    ,    0.0  , 1.0 / Vm,    0.0  },
        {     0.0,    0.0  ,      0.0    ,    0.0  ,    0.0  , 1.0 / Vm}
      }};

      mat<6,6,double> Emat = Vm * dot(transpose(Bmat), dot(inv(Jmat), Bmat));
      Emat[2][2] += Vf * Ez;

      return Emat;
    }

  };

  RotatedFiberCompositeModel(Parameters p) :
    cylindrical_domain(p.cylDom),
    stiffness_tensor(p.compute_stiffness_tensor()),
    thermal_expansion_coefficients(p.compute_thermal_expansion_coefficients()) {}

  template<typename xType, 
           typename dispType,
           typename AngleType,
           typename TemperatureType>
  auto operator()(
    const xType &x, 
    const dispType &du_dx,
    const AngleType &Alpha0,
    const TemperatureType &deltaT) {

    auto cte = thermal_expansion_coefficients;

    mat<3,3,TemperatureType> thermal_strain = {{
      {cte[0] * deltaT,             0.0,             0.0},
      {            0.0, cte[1] * deltaT,             0.0},
      {            0.0,             0.0, cte[2] * deltaT}
    }};

    // Compute rotation matrices
    using std::sin, std::cos, std::atan, std::acos, std::sqrt;
    vec3 center = {0,0,0};
    auto r = x - center; 
    auto Theta = (cylindrical_domain) * (M_PI_2 - std::atan(r[1]/r[0]));

    auto sT = sin(Theta);
    auto cT = cos(Theta);
    auto sA0 = sin(Alpha0);
    auto Alpha = acos((2.0*sA0*sT) / (sqrt(1.0 + 4.0*sA0*sA0*sT*sT))); // Dan's expression
    auto sA = sin(Alpha);
    auto cA = cos(Alpha);

    mat3 R = {{
      {cA*cT, -sT, cT*sA},
      {cA*sT,  cT, sA*sT},
      {  -sA,   0,    cA}
    }};

    vec6 wv{1.0f, 1.0f, 1.0f, 2.0f, 2.0f, 2.0f};

    auto strain = 0.5 * (du_dx + transpose(du_dx)) + thermal_strain;

    // dot(transpose(Rv), eps_v) is equivalent to:
    // 1. eps := voigt(eps_v / wv)
    // 2. eps := dot(transpose(R), eps, R)
    // 3. eps_v := wv * voigt(eps)
    strain = dot(transpose(R), dot(strain, R));

    auto stress = voigt(dot(stiffness_tensor, wv * voigt(strain)));

    // dot(Rv, sigma_v) is equivalent to:
    // 1. sigma := voigt(sigma_v)
    // 2. sigma := dot(R, sigma, transpose(R))
    // 3. sigma_v := voigt(sigma) (note: unused here)
    return dot(R, dot(stress, transpose(R)));

  }

  bool cylindrical_domain;
  mat6 stiffness_tensor; // in voigt notation
  vec3 thermal_expansion_coefficients; 

};


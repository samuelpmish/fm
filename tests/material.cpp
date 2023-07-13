#include "common.hpp"

#include <iomanip>

#include "operations/dot.hpp"
#include "operations/linear_solve.hpp"

#include "operations/random.hpp"

#include "fiber_composite/model.hpp"

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

template<typename xType, 
         typename dispType,
         typename AngleType,
         typename TemperatureType>
auto RotatedFiberComposite3D_New(
        const xType &x, 
        const dispType &du_dx,
        const AngleType &Alpha0,
        const TemperatureType &deltaT,
        double Em, double Ez, double Ep,
        double num, double nuzp, double nupp,
        double Gzp, double Vf,
        double am, double afzz, double afpp,
        const bool thermalStress, const bool cylDom) {

    // NOTE: Order used until CPlyMat is computed: 33 11 22 13 12 23 (this order is later
    // updated with a permutation matrix to match the standard Voigt notation)
    auto Vm = 1-Vf; // Matrix volume fraction

    mat<6,6,double> Jmat{};
    Jmat[0][0] = 1.0/Em;
    Jmat[0][1] = -num/Em;
    Jmat[0][2] = Jmat[0][1];
    Jmat[1][0] =  Jmat[0][1];
    Jmat[1][1] = 1.0/Em + (Vf/Vm)*(1.0/Ep - nuzp*nuzp/Ez);
    Jmat[1][2] = -num/Em - (Vf/Vm)*(nupp/Ep + nuzp*nuzp/Ez);
    Jmat[2][0] = Jmat[0][2];
    Jmat[2][1] = Jmat[1][2];
    Jmat[2][2] = Jmat[1][1];
    Jmat[3][3] = (1.0+num)/Em + (Vf/Vm)*(1.0+nupp)/Ep;
    Jmat[4][4] = (1.0+num)/Em + (Vf/Vm)*(1.0/(2.0*Gzp));
    Jmat[5][5] = Jmat[4][4];

    std::cout << std::setprecision(16);
    std::cout << "Jmat: " << Jmat << std::endl;

    auto Dmat = inv(Jmat);

    std::cout << std::setprecision(16);
    std::cout << "Dmat: " << Dmat << std::endl;

    mat<6,6,double> Bmat{};
    Bmat[0][0] = 1.0;
    Bmat[1][0] = (Vf/Vm)*nuzp;
    Bmat[2][0] = Bmat[1][0];
    Bmat[1][1] = (1.0/Vm);
    Bmat[2][2] = (1.0/Vm);
    Bmat[3][3] = (1.0/Vm);
    Bmat[4][4] = (1.0/Vm);
    Bmat[5][5] = (1.0/Vm);

    std::cout << std::setprecision(16);
    std::cout << "Bmat: " << Bmat << std::endl;

    vec<6, double> xivec{};
    xivec[1] = nuzp;
    xivec[2] = nuzp;

    mat<6, 6, double> Emat{};
    auto partialxivec = make_vec<5>([&](int i) { return xivec[i+1]; });
    auto partialDmat2 = make_vec<6>([&](int i) { return Dmat[0][i]; });
    auto partialBmat1 = make_vec<6>([&](int i) { return Bmat[i][0]; });
    auto partialDmat1 = make_mat<5, 6>([&](auto i, auto j) { return Dmat[i+1][j]; });
    auto partialBmat2 = make_mat<6, 5>([&](auto i, auto j) { return Bmat[i][j+1]; });

    // compute Emat[0][0] = Vf * (Ez + xivec(2:6) * Dmat(2:6,:) * Bmat(:,1)) + Vm * Dmat(1,:) * Bmat(:,1);
    Emat[0][0] = Vf *  (Ez + dot(partialxivec, dot(partialDmat1, partialBmat1))) + Vm * dot(partialDmat2, partialBmat1);
    // compute Emat[0][1:5] = (Vf * xivec(2:6) * Dmat(2:6,:) + Vm * Dmat(1,:)) * Bmat(:,2:6);
    auto temp1vec =  dot(Vf * dot(partialxivec, partialDmat1) + Vm * partialDmat2, partialBmat2);
    for ( int i=0; i<5; i++)
    {
        Emat[0][i+1] = temp1vec[i];
    }
    // Emat[1:5][:] = Dmat(2:6,:) * Bmat(:,:);
    auto temp1mat = dot(partialDmat1, Bmat);
    for ( int i=0; i<5; i++) {
      for ( int j=0; j<6; j++) {
        Emat[i+1][j] = temp1mat[i][j];
      }
    }

    std::cout << std::setprecision(16);
    std::cout << "Emat: " << Emat << std::endl;

    // Compute compliance
    auto Smat = inv(Emat);

    // Effective ply properties
    double G23  = 1.0 / (2.0*Smat[3][3]);
    double G13c = 1.0 / (2.0*Smat[4][4]);
    double G12c = 1.0 / (2.0*Smat[5][5]);

    double G12 = G12c;
    double G13 = G13c;

    // Correct G12 and G13 according to Hashin's estimate
    bool useHashinEstimate(false);
    if (useHashinEstimate) {
      double Gm = Em / (2.0 * (1.0 + num));
      double G12h = Gm * (Gm * Vm + Gzp * (1.0 + Vf)) / (Gm * (1.0 + Vf) + Gzp * Vm);

      // Take average of consistent and Hashin's shear moduli
      G12 = (G12c + G12h) / 2.0;
      G13 = (G13c + G12h) / 2.0;

      // Replace ply terms with correct shear terms
      // Note: no longer reordering like in LARCH:
      // Order used in this function: 11 22 33 23 31 12
      // Order that LARCH changes to (but not here): 33 11 22 12 23 31
      Smat[3][3] = 1.0/(2.0*G23); // In LARCH: 1.0/(2.0*G12);
      Smat[4][4] = 1.0/(2.0*G13); // In LARCH: 1.0/(2.0*G23);
      Smat[5][5] = 1.0/(2.0*G12); // In LARCH: 1.0/(2.0*G13);
    }

    // Compute the stiffness tensor of the ply from the compliance
    auto CPlyMat = inv(Smat); // either Smat or SPlyMat (they are the same)

    // NOTE: Despite LARCH documentation states that the order for the first
    // group of computations is zz, xx, yy, xz, xy, yz, the Matlab implementation
    // is actually xx yy zzm xy yz xz. Hence, we need to permute the stiffness tensor by
    // computing Qinv * C * S, where sigma_voigt = Q * sigma_larch; and S transforms the
    // shear components from strain to eng strain (standard voigt notation)
    mat<6,6,AngleType> SvMat{};
    SvMat[0][2] = 1.0;
    SvMat[1][0] = 1.0;
    SvMat[2][1] = 1.0;
    SvMat[3][5] = 1.0;
    SvMat[4][3] = 1.0;
    SvMat[5][4] = 1.0;

    mat<6,6,AngleType> QvMat{};
    QvMat[0][2] = 1.0;
    QvMat[1][0] = 1.0;
    QvMat[2][1] = 1.0;
    QvMat[3][5] = 0.5;
    QvMat[4][3] = 0.5;
    QvMat[5][4] = 0.5;

    auto reOrdCPlyMat = dot(transpose(SvMat), dot(CPlyMat, QvMat));

    std::cout << std::setprecision(16);
    std::cout << "reOrdCPlyMat: " << reOrdCPlyMat << std::endl;

    //// -----------------------

    // Compute rotation matrices
    using std::sin, std::cos, std::atan, std::acos, std::sqrt;
    double Xcenter = 0.0, Ycenter = 0.0;
    auto Xdist = x[0]-Xcenter;
    auto Ydist = x[1]-Ycenter;

    auto Theta = M_PI_2 - std::atan(Ydist/Xdist); // std::atan(Ydist/Xdist);
    if (!cylDom) {
      Theta = 0.0;
    }

    // For more details on variable names and operations, see the homogenize.py file of the LARCH code.
    mat<3, 3, TemperatureType > thermal_strain{};
    if (thermalStress) {

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

      auto diag = deltaT * linear_solve(Ceff, alphabar);

      thermal_strain[0][0] = diag[0];
      thermal_strain[1][1] = diag[1];
      thermal_strain[2][2] = diag[2];

    }

    vec6 wv{1.0f, 1.0f, 1.0f, 2.0f, 2.0f, 2.0f};

    auto strain = 0.5 * (du_dx + transpose(du_dx)) + thermal_strain;

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

    // dot(transpose(Rv), eps_v) is equivalent to:
    // 1. eps := voigt(eps_v / wv)
    // 2. eps := dot(transpose(R), eps, R)
    // 3. eps_v := wv * voigt(eps)
    strain = dot(transpose(R), dot(strain, R));

    auto stress = voigt(dot(reOrdCPlyMat, wv * voigt(strain)));

    // dot(Rv, sigma_v) is equivalent to:
    // 1. sigma := voigt(sigma_v)
    // 2. sigma := dot(R, sigma, transpose(R))
    // 3. sigma_v := voigt(sigma) (note: unused here)
    return dot(R, dot(stress, transpose(R)));

};

struct ModelParameters {
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
};

void run_test() {

  vec3 x = random_vec<3>();
  mat3 du_dx = random_mat<3,3>();
  double alpha0 = random_real(1.0, 3.0);
  double deltaT = random_real(1.0, 2.0);
  double Em = random_real(1.0, 3.0);
  double Ez = random_real(1.0, 3.0);
  double Ep = random_real(1.0, 3.0);
  double num = random_real(1.0, 3.0);
  double nuzp = random_real(1.0, 3.0);
  double nupp = random_real(1.0, 3.0);
  double Gzp = random_real(1.0, 3.0);
  double Vf = random_real(1.0, 3.0);
  double am = random_real(1.0, 3.0);
  double afzz = random_real(1.0, 3.0);
  double afpp = random_real(1.0, 3.0);
  bool thermalStress = true;
  bool cylDom = true;

#if 1
  std::cout << '{';
  std::cout << "deltaT, ";
  std::cout << "Em, ";
  std::cout << "Ez, ";
  std::cout << "Ep,";
  std::cout << "num,";
  std::cout << "nuzp,";
  std::cout << "nupp,";
  std::cout << "Gzp,";
  std::cout << "Vf,";
  std::cout << "am,";
  std::cout << "afzz,";
  std::cout << "afpp";
  std::cout << '}';
  std::cout << "=";
  std::cout << '{';
  std::cout << std::setprecision(16);
  std::cout << deltaT << ',';
  std::cout << Em << ',';
  std::cout << Ez << ',';
  std::cout << Ep << ',';
  std::cout << num << ',';
  std::cout << nuzp << ',';
  std::cout << nupp << ',';
  std::cout << Gzp << ',';
  std::cout << Vf << ',';
  std::cout << am << ',';
  std::cout << afzz << ',';
  std::cout << afpp;
  std::cout << '}';
  std::cout << std::endl;
#endif

  auto sigma_original = RotatedFiberComposite3D_Original(
    x, du_dx, alpha0, deltaT, 
    Em, Ez, Ep,
    num, nuzp, nupp,
    Gzp, Vf,
    am, afzz, afpp,
    thermalStress, cylDom
  );

  auto sigma_new = RotatedFiberComposite3D_New(
    x, du_dx, alpha0, deltaT, 
    Em, Ez, Ep,
    num, nuzp, nupp,
    Gzp, Vf,
    am, afzz, afpp,
    thermalStress, cylDom
  );

  compare(sigma_original, sigma_new, 5 * eps);

}

int main() {
  for (int i = 0; i < 1; i++) {
    run_test();
  }
}
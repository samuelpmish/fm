#include "common.hpp"

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
        const bool thermalStress, const bool cylDom)
{
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

    auto Dmat = inv(Jmat);

    mat<6,6,double> Bmat{};
    Bmat[0][0] = 1.0;
    Bmat[1][0] = (Vf/Vm)*nuzp;
    Bmat[2][0] = Bmat[1][0];
    Bmat[1][1] = (1.0/Vm);
    Bmat[2][2] = (1.0/Vm);
    Bmat[3][3] = (1.0/Vm);
    Bmat[4][4] = (1.0/Vm);
    Bmat[5][5] = (1.0/Vm);

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
    if (useHashinEstimate)
    {
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

    //// -----------------------

    // Compute rotation matrices
    using std::sin, std::cos,  std::atan, std::acos, std::sqrt;
    double Xcenter = 0.0, Ycenter = 0.0;
    auto Xdist = x[0]-Xcenter;
    auto Ydist = x[1]-Ycenter;

    auto Theta = M_PI_2 - std::atan(Ydist/Xdist); // std::atan(Ydist/Xdist);
    if (!cylDom) {
      Theta = 0.0;
    }

    // auto Alpha =  Alpha0;
    // auto Alpha =  M_PI_2 - acos(2.0*sin(Alpha0)*sin(Theta)/(sqrt(1.0 + 4.0*sin(Alpha0)*sin(Alpha0)*sin(Theta)*sin(Theta)))); // Dan's expression
    auto sT = sin(Theta);
    auto cT = cos(Theta);
    auto sA0 = sin(Alpha0);
    auto Alpha = acos((2.0*sA0*sT) / (sqrt(1.0 + 4.0*sA0*sA0*sT*sT))); // Dan's expression
    auto sA = sin(Alpha);
    auto cA = cos(Alpha);

    // NOTE: following pyhton LARCH implementation which is more general and performs different computations
    // than the ones above (but reached to the same final material tensor). There is some redundancy in
    // operations, but the extra cost should not be an inpediment for the size of problems solved.
    // For more details on variable names and operations, see the homogenize.py file of the LARCH code.
    mat<3, 3, TemperatureType > thermal_strain{};
    if (thermalStress)
    {
        // Generate thermal expansion coefficients
        vec<6, double> alpha_m{};
        alpha_m[0] = am;
        alpha_m[1] = am;
        alpha_m[2] = am;

        vec<6,double> alpha_f{};
        alpha_f[0] = afzz;
        alpha_f[1] = afpp;
        alpha_f[2] = afpp;

        mat<6,6,double> Smat_m{};
        Smat_m[0][0] = 1.0/Em;
        Smat_m[1][1] = 1.0/Em;
        Smat_m[2][2] = 1.0/Em;
        Smat_m[0][1] = -num/Em;
        Smat_m[0][2] = -num/Em;
        Smat_m[1][0] = -num/Em;
        Smat_m[1][2] = -num/Em;
        Smat_m[2][0] = -num/Em;
        Smat_m[2][1] = -num/Em;
        Smat_m[3][3] = (1.0+num)/Em;
        Smat_m[4][4] = (1.0+num)/Em;
        Smat_m[5][5] = (1.0+num)/Em;
        auto Cmat_m = inv(Smat_m);

        mat<6,6,double> Smat_f{};
        Smat_f[0][0] = 1.0/Ez;
        Smat_f[1][1] = 1.0/Ep;
        Smat_f[2][2] = 1.0/Ep;
        Smat_f[0][1] = -nuzp/Ez;
        Smat_f[0][2] = -nuzp/Ez;
        Smat_f[1][0] = -nuzp/Ez;
        Smat_f[2][0] = -nuzp/Ez;
        Smat_f[1][2] = -nupp/Ep;
        Smat_f[2][1] = -nupp/Ep;
        Smat_f[3][3] = 0.5/Gzp;
        Smat_f[5][5] = 0.5/Gzp;
        Smat_f[4][4] = (1.0+nupp)/Ep;
        auto Cmat_f = inv(Smat_f);

        // Extract components and perform required operations on subsets (Group 1)
        auto Cuv_m = make_mat<5, 5>([&](auto i, auto j) { return Cmat_m[i+1][j+1]; });
        auto Cuv_f = make_mat<5, 5>([&](auto i, auto j) { return Cmat_f[i+1][j+1]; });

        auto Cub_m = make_vec<5>([&](int i) { return Cmat_m[i+1][0]; });
        auto Cub_f = make_vec<5>([&](int i) { return Cmat_f[i+1][0]; });

        auto Stildn_m = inv(Cuv_m);
        auto Stildn_f = inv(Cuv_f);

        auto Btildn_m = dot(Stildn_m, Cub_m);
        auto Btildn_f = dot(Stildn_f, Cub_f);

        auto alphatildn_m = make_vec<5>([&](int i) { return alpha_m[i+1] + alpha_m[0] * Btildn_m[i]; });
        auto alphatildn_f = make_vec<5>([&](int i) { return alpha_f[i+1] + alpha_f[0] * Btildn_f[i]; });

        auto Stildbar = Vm * Stildn_m + Vf * Stildn_f;
        auto Btildbar = Vm * Btildn_m + Vf * Btildn_f;
        auto alphatildbar = make_vec<5>([&](int i) { return Vm * alphatildn_m[i] + Vf * alphatildn_f[i]; });

        auto Cbarhat = inv(Stildbar);
        auto CB = dot(Cbarhat, Btildbar);

        auto Cbarhat_Alphatildbar = dot(Cbarhat, alphatildbar);

        // Extract components and perform required operations on subsets (Group 2)
        auto Cau_m = make_vec<5>([&](int i) { return Cmat_m[0][i+1]; });
        auto Cau_f = make_vec<5>([&](int i) { return Cmat_f[0][i+1]; });

        auto Cvb_m = make_vec<5>([&](int i) { return Cmat_m[i+1][0]; });
        auto Cvb_f = make_vec<5>([&](int i) { return Cmat_f[i+1][0]; });

        auto Dtildn_m = dot(Cau_m, Stildn_m);
        auto Dtildn_f = dot(Cau_f, Stildn_f);
        auto Dtildn = Vm * Dtildn_m + Vf * Dtildn_f;

        auto Ctildn_m = dot(Dtildn_m, Cvb_m);
        auto Ctildn_f = dot(Dtildn_f, Cvb_f);

        auto Cbarbarab_m = Vm * (Cmat_m[0][0] - Ctildn_m);
        auto Cbarbarab_f = Vf * (Cmat_f[0][0] - Ctildn_f);
        auto Cbarbarab = Cbarbarab_m + Cbarbarab_f + dot(Dtildn, CB);

        auto Cbarbarav = dot(Dtildn, Cbarhat);

        auto alphabara_m = Vm * ((Cmat_m[0][0] - Ctildn_m) * alpha_m[0]);
        auto alphabara_f = Vf * ((Cmat_f[0][0] - Ctildn_f) * alpha_f[0]);
        auto alphabara = alphabara_m + alphabara_f + dot(Dtildn, alphatildbar);

        // Put all the components together into single matrices
        auto Cbarbar = make_mat<3, 3>([&](auto i, auto j)
        {
            if (i>0 && j>0) {return Cbarhat[i-1][j-1];}
            else if (i>0 && j==0) {return CB[i-1];}
            else if (i==0 && j==0) {return Cbarbarab;}
            else if (i==0 && j>0) {return Cbarbarav[j-1];}
            else {return 0.0;}
        });

        vec alphabar = {{alphabara, Cbarhat_Alphatildbar[0], Cbarhat_Alphatildbar[1]}};

        auto diag = deltaT * linear_solve(Cbarbar, alphabar);
        thermal_strain[0][0] = diag[0];
        thermal_strain[1][1] = diag[1];
        thermal_strain[2][2] = diag[2];

    }

    vec6 wv{1.0f, 1.0f, 1.0f, 2.0f, 2.0f, 2.0f};

    auto strain = 0.5 * (du_dx + transpose(du_dx)) + thermal_strain;

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
  for (int i = 0; i < 10; i++) {
    run_test();
  }
}
#pragma once

#include "fm.hpp"
#include "operations/dot.hpp"
#include "operations/linear_solve.hpp"

template<typename xType, 
         typename dispType,
         typename AngleType,
         typename TemperatureType>
auto RotatedFiberComposite3D_Original(
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
    for ( int i=0; i<5; i++)
    {
        for ( int j=0; j<6; j++)
        {
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

    bool printProperties(false);
    if (printProperties)
    {
        // Remaining effective ply properties
        double E11 = 1.0/Smat[0][0];
        double E22 = 1.0/Smat[1][1];
        double E33 = 1.0/Smat[2][2];
        double nu12 = -Smat[1][0]*E11;
        double nu13 = -Smat[2][0]*E11;
        double nu23 = -Smat[2][1]*E22;

        // Print ply properties
        std::cout<<"... E11 = "<<E11<<std::endl;
        std::cout<<"... E22 = "<<E22<<std::endl;
        std::cout<<"... E33 = "<<E33<<std::endl;
        std::cout<<"... nu12 = "<<nu12<<std::endl;
        std::cout<<"... nu13 = "<<nu13<<std::endl;
        std::cout<<"... nu23 = "<<nu23<<std::endl;
        std::cout<<"... G12 = "<<G12<<std::endl;
        std::cout<<"... G13 = "<<G13<<std::endl;
        std::cout<<"... G23 = "<<G23<<std::endl;

        std::cout<<"\n... CPlyMat = "<<inv(Smat)<<std::endl;
        std::cout<<"\n\n\n... SPlyMat = "<<Smat<<std::endl;
        exit(0);
    }

    // Compute the stiffness tensor of the ply from the compliance
    auto CPlyMat = inv(Smat); // either Smat or SPlyMat (they are the same)

    //// -----------------------
    // For debugging (using the stiffness tensor in the Appendix of Fernandez, F., Compel, W. S., Lewicki, J. P., & Tortorelli, D. A. (2019). Optimal design of fiber reinforced composite structures and their direct ink write fabrication. Computer Methods in Applied Mechanics and Engineering, 353, 277-307.)
    //// -----------------------
    bool debugFlag(false);

    // If using debug flag then the order of C is already the right voigt, i.e., xx, yy, zz, yz, xz, and xy
    // Note: engineering shear strains are used to get the stress
    if (debugFlag)
    {
        for ( int i=0; i<6; i++) {
          for ( int j=0; j<6; j++) {
            CPlyMat[i][j] = 0.0;
          }
        }

        double C11 = 9.4277;
        double C22 = 6.1028;
        double C12 = 2.9511;
        double C23 = 4.7524;
        double C44 = 1.7000;
        double C55 = 1.6300;

        CPlyMat[0][0]=C11;
        CPlyMat[0][1]=C12;
        CPlyMat[0][2]=C12;
        CPlyMat[1][0]=C12;
        CPlyMat[1][1]=C22;
        CPlyMat[1][2]=C23;
        CPlyMat[2][0]=C12;
        CPlyMat[2][1]=C23;
        CPlyMat[2][2]=C22;
        CPlyMat[3][3]=C44;
        CPlyMat[4][4]=C55;
        CPlyMat[5][5]=C44;
    }

    // NOTE: Despite LARCH documentation states that the order for the first
    // group of computations is zz, xx, yy, xz, xy, yz, the Matlab implementation
    // is actually xx yy zzm xy yz xz. Hence, we need to permutate the stiffness tensor by
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

    if (!cylDom)
    {
        Theta = 0.0;
    }

    // auto Alpha =  Alpha0;
    // auto Alpha =  M_PI_2 - acos(2.0*sin(Alpha0)*sin(Theta)/(sqrt(1.0 + 4.0*sin(Alpha0)*sin(Alpha0)*sin(Theta)*sin(Theta)))); // Dan's expression
    auto Alpha =  acos(2.0*sin(Alpha0)*sin(Theta)/(sqrt(1.0 + 4.0*sin(Alpha0)*sin(Alpha0)*sin(Theta)*sin(
                           Theta)))); // Dan's expression

    auto sA = sin(Alpha);
    auto cA = cos(Alpha);
    auto sT = sin(Theta);
    auto cT = cos(Theta);

    auto finalCPlyMat = reOrdCPlyMat;

    auto s2T = 2.0*sT*cT;
    auto c2T = cT*cT-sT*sT;
    auto s2A = 2.0*sA*cA;
    auto c2A = cA*cA-sA*sA;

    mat<6,6,AngleType> TwoRotMatVoigt{};
    TwoRotMatVoigt[0][0] = cA*cA * cT*cT;
    TwoRotMatVoigt[0][1] = sT*sT;
    TwoRotMatVoigt[0][2] = sA*sA * cT*cT;
    TwoRotMatVoigt[0][3] = -sA * s2T;
    TwoRotMatVoigt[0][4] = s2A * cT*cT;
    TwoRotMatVoigt[0][5] = -cA * s2T;

    TwoRotMatVoigt[1][0] = cA*cA * sT*sT;
    TwoRotMatVoigt[1][1] = cT*cT;
    TwoRotMatVoigt[1][2] = sA*sA * sT*sT;
    TwoRotMatVoigt[1][3] = sA * s2T;
    TwoRotMatVoigt[1][4] = s2A * sT*sT;
    TwoRotMatVoigt[1][5] = cA * s2T;

    TwoRotMatVoigt[2][0] = sA*sA;
    TwoRotMatVoigt[2][1] = 0.0;
    TwoRotMatVoigt[2][2] = cA*cA;
    TwoRotMatVoigt[2][3] = 0.0;
    TwoRotMatVoigt[2][4] = -s2A;
    TwoRotMatVoigt[2][5] = 0.0;

    TwoRotMatVoigt[3][0] = -0.5*s2A * sT;
    TwoRotMatVoigt[3][1] = 0.0;
    TwoRotMatVoigt[3][2] = 0.5*s2A * sT;
    TwoRotMatVoigt[3][3] = cA* cT;
    TwoRotMatVoigt[3][4] = c2A * sT;
    TwoRotMatVoigt[3][5] = -sA * cT;

    TwoRotMatVoigt[4][0] = -0.5*s2A * cT;
    TwoRotMatVoigt[4][1] = 0.0;
    TwoRotMatVoigt[4][2] = 0.5*s2A * cT;
    TwoRotMatVoigt[4][3] = -cA * sT;
    TwoRotMatVoigt[4][4] = c2A * cT;
    TwoRotMatVoigt[4][5] = sA * sT;

    TwoRotMatVoigt[5][0] = 0.5*cA*cA * s2T;
    TwoRotMatVoigt[5][1] = -0.5*s2T;
    TwoRotMatVoigt[5][2] = 0.5*sA*sA * s2T;
    TwoRotMatVoigt[5][3] = sA * c2T;
    TwoRotMatVoigt[5][4] = 0.5*s2A * s2T;
    TwoRotMatVoigt[5][5] = cA * c2T;

    // include effect of rotation
    finalCPlyMat = dot(TwoRotMatVoigt, dot(reOrdCPlyMat, transpose(TwoRotMatVoigt)));

    // NOTE: following pyhton LARCH implementation which is more general and performs different computations
    // than the ones above (but reached to the same final material tensor). There is some redundancy in
    // operations, but the extra cost should not be an inpediment for the size of problems solved.
    // For more details on variable names and operations, see the homogenize.py file of the LARCH code.
    vec<6, TemperatureType > thermalStrainVec{};
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

        auto Ctildn_m = dot(Dtildn_m, Cvb_m);
        auto Ctildn_f = dot(Dtildn_f, Cvb_f);

        auto Cbarbarab_m = Vm * (Cmat_m[0][0] - Ctildn_m + dot(Dtildn_m, CB));
        auto Cbarbarab_f = Vf * (Cmat_f[0][0] - Ctildn_f + dot(Dtildn_f, CB));
        auto Cbarbarab = Cbarbarab_m + Cbarbarab_f;

        auto Cbarbarav_m = Vm * dot(Dtildn_m, Cbarhat);
        auto Cbarbarav_f = Vf * dot(Dtildn_f, Cbarhat);
        auto Cbarbarav = Cbarbarav_m + Cbarbarav_f;

        auto alphabara_m = Vm * ((Cmat_m[0][0]-Ctildn_m) *  alpha_m[0] + dot(Dtildn_m, alphatildbar));
        auto alphabara_f = Vf * ((Cmat_f[0][0]-Ctildn_f) *  alpha_f[0] + dot(Dtildn_f, alphatildbar));
        auto alphabara = alphabara_m + alphabara_f;

        // Put all the components together into single matrices
        auto Cbarbar = make_mat<6, 6>([&](auto i, auto j)
        {
            if (i>0 && j>0) {return Cbarhat[i-1][j-1];}
            else if (i>0 && j==0) {return CB[i-1];}
            else if (i==0 && j==0) {return Cbarbarab;}
            else if (i==0 && j>0) {return Cbarbarav[j-1];}
            else {return 0.0;}
        });

        auto alphabar = make_vec<6>([&](int i)
        {
            if (i==0) {return alphabara;}
            else if (i>0) {return Cbarhat_Alphatildbar[i-1];}
            else {return 0.0;}
        });

        // Compute the thermal expansion coefficient for the ply
        auto Sbarbar = inv(Cbarbar);
        auto finalAlphaPlyVec = dot(Sbarbar, alphabar);

        // Compute thermal stress
        thermalStrainVec = deltaT * finalAlphaPlyVec;

    }

    // Compute stress
    auto strain = 0.5 * (du_dx + transpose(du_dx)); // small strain tensor

    // vectorize tensor using the order: xx, yy, zz, 2yz, 2xz, 2xy
    auto strucStrainVec = make_vec<6>([&](int i)
    {
        if (i<3) {return strain(i,i);}
        else if (i==3) {return 2.0*strain(1,2);}
        else if (i==4) {return 2.0*strain(0,2);}
        else {return 2.0*strain(0,1);}
    });

    auto totalStrainVec = strucStrainVec + thermalStrainVec;

    auto totalStressVec = dot(finalCPlyMat, totalStrainVec);

    // if(totalStressVec[0]>0.0)
    // {
    // std::cout<<".................."<<std::endl;
    // std::cout<< serac::get_value(totalStressVec) << std::endl;
    // std::cout<<".................."<<std::endl;
    // std::cout<< thermalStrainVec << std::endl;
    // std::cout<<".................."<<std::endl;
    // std::cout<< serac::get_value(totalStrainVec) << std::endl;
    // std::cout<<".................."<<std::endl;
    // exit(0);
    // }

    // matricize stress because return expects a 3x3 tensor
    auto stress = make_mat<3,3>([&](auto i, auto j)
    {
        if (i==j) {return totalStressVec(i);}
        else if ((i==1&&j==2) || (i==2&&j==1)) {return totalStressVec(3);}
        else if ((i==0&&j==2) || (i==2&&j==0)) {return totalStressVec(4);}
        else {return totalStressVec(5);}
    });

    return stress;
};
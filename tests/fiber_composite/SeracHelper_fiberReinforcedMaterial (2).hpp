// Copyright (c) 2022, Lawrence Livermore National Security, LLC.
// All rights reserved.  LLNL-CODE-728517

// OFFICIAL USE ONLY This work was produced at the Lawrence Livermore
// National Laboratory (LLNL) under contract no. DE-AC52-07NA27344
// (Contract 44) between the U.S. Department of Energy (DOE) and
// Lawrence Livermore National Security, LLC (LLNS) for the operation
// of LLNL.  See license for disclaimers, notice of U.S. Government
// Rights and license terms and conditions.

#include "lido_config.hpp"

#ifndef LIDO_USE_SERAC
static_assert(false, "LiDO was not built with serac support");
#else // defined(LIDO_USE_SERAC)

#ifndef __LIDO_SERAC_HELPER_HPP__
#define __LIDO_SERAC_HELPER_HPP__

#include <mfem.hpp>
#include <serac/numerics/functional/functional.hpp>
#include <serac/physics/state/finite_element_state.hpp>

#define TWO_ROT_IN_ONE
// #undef TWO_ROT_IN_ONE

namespace lido
{

/**
 * @brief A factory function to generate a serac::FiniteElementState on a templated finite element space
 *
 * @tparam space    Finite element space, e.g. serac::H1<p, vdim> for an H1 space of order p with vector dimension vdim
 *                  A quadratic scalar space would be serac::H1<2, 1>, whereas a linear vector space could be serac::H1<1, 3>
 * @param pmesh     Reference to parallel mesh
 * @param name      Description of space (must be unique amongst FiniteElementStates registered with serac)
 * @return          std::unique_ptr<serac::FiniteElementState>
 */
template<typename space>
std::unique_ptr<serac::FiniteElementState> StateFactory(mfem::ParMesh &pmesh, std::string name)
{
    const int poly_order = space::order;
    const int vector_dim = space::components;

    if (space::family == serac::Family::L2) // return an L2 finite element state
    {
        return std::make_unique<serac::FiniteElementState>(pmesh, serac::FiniteElementState::Options
        {
            .order = poly_order,
            .vector_dim=vector_dim,
            .element_type=serac::ElementType::L2,
            .name=name});
    }
    else if (space::family == serac::Family::H1) // return an H1 finite element state
    {
        return std::make_unique<serac::FiniteElementState>(pmesh, serac::FiniteElementState::Options
        {
            .order = poly_order,
            .vector_dim=vector_dim,
            .element_type=serac::ElementType::H1,
            .name=name});
    }
    else // throw if T::family is anything other than L2 or H1, e.g. HCurl or HDiv
    {
        throw GraphArgumentException("Only L2 and H1 FE spaces are supported for serac parameters", __FILE__, __LINE__);
    }
}

/**
 * @brief   A material model for use in serac::solidFunctional.
 *          This model is a linear, isotropic material defined in terms of
 *          Young's modulus and Poisson's ratio
 * @details Materials in serac always require 3 things-
 *          1. A material @p State
 *          2. An @p () overload that computes the stress
 *          3. The @p density (only used in dynamics problems)
 */
struct LinearIsotropicSolid
{
    using State = serac::Empty;  ///< this material has no internal variables

    /**
     * @brief Construct a new LinearIsotropicSolid_Lame material
     *
     * @param E Young's modulus
     * @param v Poisson's ratio
     */
    LinearIsotropicSolid(double E, double v)
        : lambda_(E*v / ((1.0+v)*(1.0-2.0*v))),
          mu_(E / (2.0*(1.0+v)))
    {}

    /// default constructor not allowed
    LinearIsotropicSolid() = delete;

    /**
     * @brief  The operator overload must return stress as a function of material state,
     *         displacement gradient, and any material parameters
     *
     * @tparam dim           The spatial dimension of the problem
     * @tparam T             The type of the displacement gradient (serac uses this internally)
     * @param u_grad         The displacement spatial gradient
     */
    template <int dim, typename T>
    SERAC_HOST_DEVICE auto operator()(State &, const serac::tensor<T, dim, dim> &u_grad) const
    {
        const auto I      = serac::Identity<dim>();       // identity matrix
        auto       strain = serac::sym(u_grad);           // small strain tensor
        return lambda_*serac::tr(strain)*I + 2.0*mu_*strain;  // Cauchy stress
    }

    static constexpr double density{1.0}; ///< mass density, for dynamics problems

    double lambda_, mu_;
};

/**
 * @brief   A parameterized material model for use in serac::solidFunctional.
 *          This model is a linear, isotropic material defined in terms of
 *          Young's modulus and Poisson's ratio.
 * @details Materials in serac always require 3 things-
 *          1. A material @p State
 *          2. An @p () overload that computes the stress
 *          3. The @p density (only used in dynamics problems)
 */
struct ParameterizedLinearIsotropicSolid
{
    using State = serac::Empty;  ///< this material has no internal variables

    /**
     * @brief  The operator overload must return stress as a function of material state,
     *         displacement gradient, and any material parameters
     *
     * @tparam dim           The spatial dimension of the problem
     * @tparam T             The type of the displacement gradient (serac uses this internally)
     * @tparam YoungsType    The type of the Young's modulus
     * @tparam PoissonsType  The type of the Poisson's ratio
     * @param u_grad         The displacement spatial gradient
     * @param E_tuple        Young's modulus and its gradient
     * @param v_tuple        Poisson's ratio and its gradient
     */
    template <int dim, typename T, typename YoungsType, typename PoissonsType>
    SERAC_HOST_DEVICE auto operator()(State &, const serac::tensor<T, dim, dim> &u_grad,
                                      const YoungsType &E_tuple, const PoissonsType &v_tuple) const
    {
        auto       E      = serac::get<0>(E_tuple);      // Young's modulus VALUE
        auto       v      = serac::get<0>(v_tuple);      // Poisson's ratio VALUE
        auto       lambda = E*v / ((1.0+v)*(1.0-2.0*v)); // Lamé's first parameter
        auto       mu     = E / (2.0*(1.0+v));           // Lamé's second parameter
        const auto I      = serac::Identity<dim>();      // identity matrix
        auto       strain = serac::sym(u_grad);          // small strain tensor
        return lambda*serac::tr(strain)*I + 2.0*mu*strain;   // Cauchy stress
    }

    static constexpr double density{1.0}; ///< mass density, for dynamics problems
};

/**
 * @brief A compliance QoI that assumes constant material properties, but variable shape
 */
struct ShapeAwareLinearIsotropicCompliance
{
    /**
     * @brief Construct a new ShapeAwareLinearIsotropicCompliance object
     *
     * @param E  Young's modulus
     * @param v  Poisson's ratio
     */
    ShapeAwareLinearIsotropicCompliance(double E, double v)
        :  material_(E, v)
    {}

    /// default constructor not allowed
    ShapeAwareLinearIsotropicCompliance() = delete;

    /**
     * @brief The operator overload returns compliance as a function of shape and displacement
     *
     * @tparam dim           The spatial dimension of the problem
     * @tparam T             The type of the displacement gradient (serac uses this internally)
     * @tparam shapeType     The type of the shape field
     * @tparam dispType      The type of the displacement field
     * @param x              Spatial location
     * @param X_tuple        The shape field and its gradient
     * @param u_tuple        The displacement and its gradient
     */
    template <int dim, typename T, typename shapeType, typename dispType>
    SERAC_HOST_DEVICE auto operator()(const serac::tensor<T, dim> LIDO_UNUSED(x), const shapeType &X_tuple,
                                      const dispType &u_tuple)
    {
        auto       state     = serac::Empty{};                         // assuming state-free material
        auto       dp_dx     = serac::get<1>(X_tuple);                 // gradient of shape field
        const auto I         = serac::Identity<dim>();                 // identity tensor
        auto       du_dx     = serac::get<1>
                               (u_tuple);                 // compute displacement gradient (w.r.t. perturbed domain)
        auto       du_dx_ref = serac::dot(du_dx,
                                          serac::inv(I+dp_dx)); // compute displacement gradient (w.r.t. reference domain)
        auto       strain    = serac::sym(du_dx_ref);                  // small strain tensor
        auto       stress    = material_(state, du_dx_ref);            // cauchy stress
        return 0.5 * serac::double_dot(strain, stress) * serac::det(dp_dx + I);
    }

    LinearIsotropicSolid material_;
};

template<typename xType, typename dispType, typename YoungsType, typename PoissonsType, typename ShearType>
SERAC_HOST_DEVICE auto RotatedCubicStress3D(const xType &x, const dispType &u, const YoungsType &E_tuple,
        const PoissonsType &v_tuple, const ShearType &G_tuple)
{
    using std::sin, std::cos, std::atan2;

    auto E     = serac::get<0>(E_tuple);  // Young's modulus VALUE
    auto v     = serac::get<0>(v_tuple);  // Poisson's ratio VALUE
    auto G     = serac::get<0>(G_tuple);  // shear modulus VALUE
    auto c1111 = E*(1-v) / (1-v-2*v*v);   // cubic constant
    auto c1122 = E*v / ((1+v)*(1-2*v));   // cubic constant (Lamé's first parameter)
    auto c1212 = G;                       // cubic constant (Lamé's second parameter)

    auto y = x[1] - 0.179528;
    auto z = x[2] - -5.5090945;
    auto t = atan2(z, y);

    auto c2222 = (3*c1111 +   c1122 + 2*c1212 + cos(4*t)*( c1111 - c1122 - 2*c1212)) / 4;
    auto c2233 = (  c1111 + 3*c1122 - 2*c1212 + cos(4*t)*(-c1111 + c1122 + 2*c1212)) / 4;
    auto c3232 = (  c1111 -   c1122 + 2*c1212 + cos(4*t)*(-c1111 + c1122 + 2*c1212)) / 4;
    auto c2232 =  (sin(4*t)*(c1111 - c1122 - 2*c1212)) / 4;

    auto du_dx  = serac::get<1>(u);     // grab displacement GRADIENT
    auto strain = serac::sym(du_dx);    // small strain tensor
    auto stress = 0.0 * strain * E * v * G;

    stress[0][0] += c1111*strain[0][0] + c1122*strain[1][1] + c1122*strain[2][2];
    stress[1][1] += c1122*strain[0][0] + c2222*strain[1][1] + 2*c2232*strain[1][2] + c2233*strain[2][2];
    stress[2][2] += c1122*strain[0][0] + c2233*strain[1][1] - 2*c2232*strain[1][2] + c2222*strain[2][2];
    stress[0][1] += 2*c1212*strain[0][1];
    stress[0][2] += 2*c1212*strain[0][2];
    stress[1][2] += c2232*strain[1][1] + 2*c3232*strain[1][2] - c2232*strain[2][2];
    stress[1][0] += stress[0][1];
    stress[2][0] += stress[0][2];
    stress[2][1] += stress[1][2];

    return stress;
};

struct RotatedLinearCubicSolid3D
{
    template <typename xType, typename dispType, typename accelType, typename shapeType, typename YoungsType, typename PoissonsType, typename ShearType>
    SERAC_HOST_DEVICE auto operator()(const xType &x, const dispType &u, const accelType &a, const shapeType &,
                                      const YoungsType &E_tuple, const PoissonsType &v_tuple, const ShearType &G_tuple) const
    {
        auto d2u_dt2 = serac::get<0>(a);
        auto stress  = RotatedCubicStress3D(x, u, E_tuple, v_tuple, G_tuple);
        return serac::tuple{d2u_dt2, stress};
    }
};

template<typename xType, typename dispType,
         typename AngleType, typename AngleGradType,
         typename TemperatureType, typename TemperatureGradType>
SERAC_HOST_DEVICE auto RotatedFiberComposite3D(const xType &x, const dispType &u,
        const serac::tuple< AngleType, AngleGradType > &Alpha_tuple,
        const serac::tuple< TemperatureType, TemperatureGradType > &T,
        double Em, double Ez, double Ep,
        double num, double nuzp, double nupp,
        double Gzp, double Vf,
        double am, double afzz, double afpp,
        const bool thermalStress, const bool cylDom)
{
    // NOTE: Order used until CPlyMat is computed: 33 11 22 13 12 23 (this order is later
    // updated with a permutation matrix to match the standard Voigt notation)
    auto Vm = 1-Vf; // Matrix volume fraction

    serac::tensor<double, 6, 6> Jmat{};
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

    serac::tensor<double, 6, 6> Bmat{};
    Bmat[0][0] = 1.0;
    Bmat[1][0] = (Vf/Vm)*nuzp;
    Bmat[2][0] = Bmat[1][0];
    Bmat[1][1] = (1.0/Vm);
    Bmat[2][2] = (1.0/Vm);
    Bmat[3][3] = (1.0/Vm);
    Bmat[4][4] = (1.0/Vm);
    Bmat[5][5] = (1.0/Vm);

    serac::tensor<double, 6> xivec{};
    xivec[1] = nuzp;
    xivec[2] = nuzp;

    serac::tensor<double, 6, 6> Emat{};
    auto partialxivec = serac::make_tensor<5>([&](auto i) { return xivec(i+1); });
    auto partialDmat2 = serac::make_tensor<6>([&](auto i) { return Dmat(0, i); });
    auto partialBmat1 = serac::make_tensor<6>([&](auto i) { return Bmat(i, 0); });
    auto partialDmat1 = serac::make_tensor<5, 6>([&](auto i, auto j) { return Dmat(i+1, j); });
    auto partialBmat2 = serac::make_tensor<6, 5>([&](auto i, auto j) { return Bmat(i, j+1); });

    // compute Emat[0][0] = Vf * (Ez + xivec(2:6) * Dmat(2:6,:) * Bmat(:,1)) + Vm * Dmat(1,:) * Bmat(:,1);
    Emat[0][0] = Vf *  (Ez + dot(partialxivec, dot(partialDmat1, partialBmat1))) + Vm * dot(partialDmat2, partialBmat1);
    // compute Emat[0][1:5] = (Vf * xivec(2:6) * Dmat(2:6,:) + Vm * Dmat(1,:)) * Bmat(:,2:6);
    auto temp1vec =  dot(Vf * dot(partialxivec, partialDmat1) + Vm * partialDmat2, partialBmat2);
    for ( int i=0; i<5; i++)
    {
        Emat[0][i+1] = temp1vec(i);
    }
    // Emat[1:5][:] = Dmat(2:6,:) * Bmat(:,:);
    auto temp1mat = dot(partialDmat1, Bmat);
    for ( int i=0; i<5; i++)
    {
        for ( int j=0; j<6; j++)
        {
            Emat[i+1][j] = temp1mat(i,j);
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
        for ( int i=0; i<6; i++)
            for ( int j=0; j<6; j++)
            {
                CPlyMat(i,j) = 0.0;
            }

        double C11 = 9.4277;
        double C22 = 6.1028;
        double C12 = 2.9511;
        double C23 = 4.7524;
        double C44 = 1.7000;
        double C55 = 1.6300;

        CPlyMat(0,0)=C11;
        CPlyMat(0,1)=C12;
        CPlyMat(0,2)=C12;
        CPlyMat(1,0)=C12;
        CPlyMat(1,1)=C22;
        CPlyMat(1,2)=C23;
        CPlyMat(2,0)=C12;
        CPlyMat(2,1)=C23;
        CPlyMat(2,2)=C22;
        CPlyMat(3,3)=C44;
        CPlyMat(4,4)=C55;
        CPlyMat(5,5)=C44;
    }

    // NOTE: Despite LARCH documentation states that the order for the first
    // group of computations is zz, xx, yy, xz, xy, yz, the Matlab implementation
    // is actually xx yy zzm xy yz xz. Hence, we need to permutate the stiffness tensor by
    // computing Qinv * C * S, where sigma_voigt = Q * sigma_larch; and S transforms the
    // shear components from strain to eng strain (standard voigt notation)
    serac::tensor<AngleType, 6, 6> SvMat{};
    SvMat[0][2] = 1.0;
    SvMat[1][0] = 1.0;
    SvMat[2][1] = 1.0;
    SvMat[3][5] = 1.0;
    SvMat[4][3] = 1.0;
    SvMat[5][4] = 1.0;

    serac::tensor<AngleType, 6, 6> QvMat{};
    QvMat[0][2] = 1.0;
    QvMat[1][0] = 1.0;
    QvMat[2][1] = 1.0;
    QvMat[3][5] = 0.5;
    QvMat[4][3] = 0.5;
    QvMat[5][4] = 0.5;

    auto reOrdCPlyMat = serac::transpose(SvMat) * CPlyMat * QvMat;

    //// -----------------------

    // Compute rotation matrices
    using std::sin, std::cos,  std::atan, std::acos, std::sqrt;
    double Xcenter = 0.0, Ycenter = 0.0;
    auto Xdist = x[0]-Xcenter;
    auto Ydist = x[1]-Ycenter;

    auto Alpha0 = serac::get<0>(Alpha_tuple);
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

    serac::tensor<AngleType, 6, 6> TwoRotMatVoigt{};
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
    finalCPlyMat = TwoRotMatVoigt * reOrdCPlyMat * serac::transpose(TwoRotMatVoigt);

    // NOTE: following pyhton LARCH implementation which is more general and performs different computations
    // than the ones above (but reached to the same final material tensor). There is some redundancy in
    // operations, but the extra cost should not be an inpediment for the size of problems solved.
    // For more details on variable names and operations, see the homogenize.py file of the LARCH code.
    serac::tensor< TemperatureType, 6 > thermalStrainVec{};
    if (thermalStress)
    {
        // Generate thermal expansion coefficients
        serac::tensor<double, 6> alpha_m{};
        alpha_m[0] = am;
        alpha_m[1] = am;
        alpha_m[2] = am;

        serac::tensor<double, 6> alpha_f{};
        alpha_f[0] = afzz;
        alpha_f[1] = afpp;
        alpha_f[2] = afpp;

        serac::tensor<double, 6, 6> Smat_m{};
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
        auto Cmat_m = serac::inv(Smat_m);

        serac::tensor<double, 6, 6> Smat_f{};
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
        auto Cmat_f = serac::inv(Smat_f);

        // Extract components and perform required operations on subsets (Group 1)
        auto Cuv_m = serac::make_tensor<5, 5>([&](auto i, auto j) { return Cmat_m(i+1, j+1); });
        auto Cuv_f = serac::make_tensor<5, 5>([&](auto i, auto j) { return Cmat_f(i+1, j+1); });

        auto Cub_m = serac::make_tensor<5>([&](auto i) { return Cmat_m(i+1, 0); });
        auto Cub_f = serac::make_tensor<5>([&](auto i) { return Cmat_f(i+1, 0); });

        auto Stildn_m = serac::inv(Cuv_m);
        auto Stildn_f = serac::inv(Cuv_f);

        auto Btildn_m = Stildn_m * Cub_m;
        auto Btildn_f = Stildn_f * Cub_f;

        auto alphatildn_m = serac::make_tensor<5>([&](auto i) { return alpha_m(i+1) + alpha_m[0] * Btildn_m(i); });
        auto alphatildn_f = serac::make_tensor<5>([&](auto i) { return alpha_f(i+1) + alpha_f[0] * Btildn_f(i); });

        auto Stildbar = Vm * Stildn_m + Vf * Stildn_f;
        auto Btildbar = Vm * Btildn_m + Vf * Btildn_f;
        auto alphatildbar = serac::make_tensor<5>([&](auto i) { return Vm * alphatildn_m(i) + Vf * alphatildn_f(i); });

        auto Cbarhat = serac::inv(Stildbar);
        auto CB = Cbarhat * Btildbar;

        auto Cbarhat_Alphatildbar = Cbarhat * alphatildbar;

        // Extract components and perform required operations on subsets (Group 2)
        auto Cau_m = serac::make_tensor<5>([&](auto i) { return Cmat_m(0, i+1); });
        auto Cau_f = serac::make_tensor<5>([&](auto i) { return Cmat_f(0, i+1); });

        auto Cvb_m = serac::make_tensor<5>([&](auto i) { return Cmat_m(i+1, 0); });
        auto Cvb_f = serac::make_tensor<5>([&](auto i) { return Cmat_f(i+1, 0); });

        auto Dtildn_m = Cau_m * Stildn_m;
        auto Dtildn_f = Cau_f * Stildn_f;

        auto Ctildn_m = serac::dot(Dtildn_m, Cvb_m);
        auto Ctildn_f = serac::dot(Dtildn_f, Cvb_f);

        auto Cbarbarab_m = Vm * (Cmat_m[0][0] - Ctildn_m + serac::dot(Dtildn_m, CB));
        auto Cbarbarab_f = Vf * (Cmat_f[0][0] - Ctildn_f + serac::dot(Dtildn_f, CB));
        auto Cbarbarab = Cbarbarab_m + Cbarbarab_f;

        auto Cbarbarav_m = Vm * serac::dot(Dtildn_m, Cbarhat);
        auto Cbarbarav_f = Vf * serac::dot(Dtildn_f, Cbarhat);
        auto Cbarbarav = Cbarbarav_m + Cbarbarav_f;

        auto alphabara_m = Vm * ((Cmat_m[0][0]-Ctildn_m) *  alpha_m[0] + serac::dot(Dtildn_m, alphatildbar));
        auto alphabara_f = Vf * ((Cmat_f[0][0]-Ctildn_f) *  alpha_f[0] + serac::dot(Dtildn_f, alphatildbar));
        auto alphabara = alphabara_m + alphabara_f;

        // Put all the components together into single matrices
        auto Cbarbar = serac::make_tensor<6, 6>([&](auto i, auto j)
        {
            if (i>0 && j>0) {return Cbarhat(i-1, j-1);}
            else if (i>0 && j==0) {return CB(i-1);}
            else if (i==0 && j==0) {return Cbarbarab;}
            else if (i==0 && j>0) {return Cbarbarav(j-1);}
            else {return 0.0;}
        });

        auto alphabar = serac::make_tensor<6>([&](auto i)
        {
            if (i==0) {return alphabara;}
            else if (i>0) {return Cbarhat_Alphatildbar(i-1);}
            else {return 0.0;}
        });

        // Compute the thermal expansion coefficient for the ply
        auto Sbarbar = serac::inv(Cbarbar);
        auto finalAlphaPlyVec = Sbarbar * alphabar;

        // Compute thermal stress
        auto deltaT  = serac::get<0>(T);
        thermalStrainVec = deltaT * finalAlphaPlyVec;

    }

    // Compute stress
    auto du_dx  = serac::get<1>(u);     // grab displacement GRADIENT
    auto strain = serac::sym(du_dx);    // small strain tensor

    // vectorize tensor using the order: xx, yy, zz, 2yz, 2xz, 2xy
    auto strucStrainVec = serac::make_tensor<6>([&](auto i)
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
    auto stress = serac::make_tensor<3,3>([&](auto i, auto j)
    {
        if (i==j) {return totalStressVec(i);}
        else if ((i==1&&j==2) || (i==2&&j==1)) {return totalStressVec(3);}
        else if ((i==0&&j==2) || (i==2&&j==0)) {return totalStressVec(4);}
        else {return totalStressVec(5);}
    });

    return stress;
};

struct RotatedLinearFiberComposite3D
{
    RotatedLinearFiberComposite3D(double Em, double Ez, double Ep,
                                  double num, double nuzp, double nupp, double Gzp, double Vf,
                                  double am, double afzz, double afpp,
                                  bool thermalStress=false, bool cylDom=true):
        Em_(Em), Ez_(Ez), Ep_(Ep), num_(num), nuzp_(nuzp), nupp_(nupp), Gzp_(Gzp), Vf_(Vf),
        am_(am), afzz_(afzz), afpp_(afpp), thermalStress_(thermalStress), cylDom_(cylDom) {}

    RotatedLinearFiberComposite3D() = delete;

    template <typename xType, typename dispType, typename accelType, typename shapeType,
              typename AngleType, typename AngleGradType,
              typename TemperatureType, typename TemperatureGradType>
    SERAC_HOST_DEVICE auto operator()(
        const xType &x, const dispType &u, const accelType &a, const shapeType &,
        const serac::tuple< AngleType, AngleGradType > &Alpha_tuple,
        const serac::tuple< TemperatureType, TemperatureGradType > deltaTemp) const
    {
        auto d2u_dt2 = serac::get<0>(a);
        auto stress  = RotatedFiberComposite3D(x, u, Alpha_tuple, deltaTemp,
                                               Em_, Ez_, Ep_, num_, nuzp_, nupp_, Gzp_, Vf_,
                                               am_, afzz_, afpp_, thermalStress_, cylDom_);
        return serac::tuple{d2u_dt2, stress};
    }
    double Em_, Ez_, Ep_, num_, nuzp_, nupp_, Gzp_, Vf_, am_, afzz_, afpp_;
    bool thermalStress_, cylDom_;
};

/**
 * @brief A compliance QoI that assumes constant shape, but variable material properties
 */
struct MaterialAwareLinearIsotropicCompliance
{
    /**
     * @brief The operator overload returns compliance as a function of material properties and displacement
     *
     * @tparam dim           The spatial dimension of the problem
     * @tparam T             The type of the displacement gradient (serac uses this internally)
     * @tparam YoungsType    The type of the Young's modulus
     * @tparam PoissonsType  The type of the Poissons' ratio
     * @tparam dispType      The type of the displacement
     * @param x              Spatial location
     * @param E_tuple        Young's modulus and its gradient
     * @param v_tuple        Poisson's ratio and its gradient
     * @param u_tuple        Displacement and its gradient
     */
    template <int dim, typename T, typename YoungsType, typename PoissonsType, typename dispType>
    SERAC_HOST_DEVICE auto operator()(const serac::tensor<T, dim> LIDO_UNUSED(x), const YoungsType &E_tuple,
                                      const PoissonsType &v_tuple, const dispType &u_tuple)
    {
        auto state  = serac::Empty{};          // assuming state-free material
        auto du_dx  = serac::get<1>(u_tuple);  // compute displacement gradient
        auto strain = serac::sym(du_dx);       // small strain tensor
        auto stress = material_(state, du_dx, E_tuple, v_tuple); // cauchy stress
        return 0.5 * serac::double_dot(strain, stress);
    }

    ParameterizedLinearIsotropicSolid material_;
};

/**
 * @brief A compliance QoI that assumes variable shape and material properties
 */
struct ShapeAndMaterialAwareLinearIsotropicCompliance
{
    /**
     * @brief The operator overload returns compliance as a function of shape, material properties, and displacement
     *
     * @tparam dim           The spatial dimension of the problem
     * @tparam T             The type of the displacement gradient (serac uses this internally)
     * @tparam shapeType     The type of the shape field
     * @tparam YoungsType    The type of the Young's modulus
     * @tparam PoissonsType  The type of the Poissons' ratio
     * @tparam dispType      The type of the displacement
     * @param x              Spatial location
     * @param X_tuple        The shape field and its gradient
     * @param E_tuple        Young's modulus and its gradient
     * @param v_tuple        Poisson's ratio and its gradient
     * @param u_tuple        Displacement and its gradient
     */
    template <int dim, typename T, typename shapeType, typename YoungsType, typename PoissonsType, typename dispType>
    SERAC_HOST_DEVICE auto operator()(const serac::tensor<T, dim> LIDO_UNUSED(x), const shapeType &X_tuple,
                                      const YoungsType &E_tuple, const PoissonsType &v_tuple, const dispType &u_tuple)
    {
        auto       state     = serac::Empty{};                         // assuming state-free material
        auto       dp_dx     = serac::get<1>(X_tuple);                 // gradient of shape field
        const auto I         = serac::Identity<dim>();                 // identity tensor
        auto       du_dx     = serac::get<1>
                               (u_tuple);                 // compute displacement gradient (w.r.t. perturbed domain)
        auto       du_dx_ref = serac::dot(du_dx,
                                          serac::inv(I+dp_dx)); // compute displacement gradient (w.r.t. reference domain)
        auto       strain    = serac::sym(du_dx_ref);                  // small strain tensor
        auto       stress    = material_(state, du_dx_ref, E_tuple, v_tuple);  // cauchy stress
        return 0.5 * serac::double_dot(strain, stress) * serac::det(dp_dx + I);
    }

    ParameterizedLinearIsotropicSolid material_;
};

/**
 * @brief A mass QoI that assumes constant density (1.0)
 */
struct ShapeAwareMass
{
    /**
     * @brief The operator overload returns shape-adjusted density
     *
     * @tparam dim           The spatial dimension of the problem
     * @tparam T             The type of the displacement gradient (serac uses this internally)
     * @tparam shapeType     The type of the shape field
     * @param x              Spatial location
     * @param X_tuple        The shape field and its gradient
     */
    template <int dim, typename T, typename shapeType>
    SERAC_HOST_DEVICE auto operator()(const serac::tensor<T, dim> LIDO_UNUSED(x), const shapeType &X_tuple)
    {
        auto       dp_dx = serac::get<1>(X_tuple);  // gradient of shape field
        const auto I     = serac::Identity<dim>();  // identity tensor
        return serac::det(dp_dx + I);
    }
};

/**
 * @brief A mass QoI that assumes variable density (but constant shape)
 */
struct MaterialAwareMass
{
    /**
     * @brief The operator overload returns density
     *
     * @tparam dim           The spatial dimension of the problem
     * @tparam T             The type of the displacement gradient (serac uses this internally)
     * @tparam densType      The type of the density
     * @param x              Spatial location
     * @param rho_tuple      The density and its gradient
     */
    template <int dim, typename T, typename densType>
    SERAC_HOST_DEVICE auto operator()(const serac::tensor<T, dim> LIDO_UNUSED(x), const densType &rho_tuple)
    {
        return serac::get<0>(rho_tuple);
    }
};

/**
 * @brief A mass QoI that assumes variable shape and density
 */
struct ShapeAndMaterialAwareMass
{
    /**
     * @brief The operator overload returns shape-adjusted, variable density
     *
     * @tparam dim           The spatial dimension of the problem
     * @tparam T             The type of the displacement gradient (serac uses this internally)
     * @tparam shapeType     The type of the shape field
     * @tparam densType      The type of the density
     * @param x              Spatial location
     * @param X_tuple        The shape field and its gradient
     * @param rho_tuple      The density and its gradient
     */
    template <int dim, typename T, typename shapeType, typename densType>
    SERAC_HOST_DEVICE auto operator()(const serac::tensor<T, dim> LIDO_UNUSED(x), const shapeType &X_tuple,
                                      const densType &rho_tuple)
    {
        auto       dp_dx = serac::get<1>(X_tuple);  // gradient of shape field
        const auto I     = serac::Identity<dim>();  // identity tensor
        return serac::get<0>(rho_tuple) * serac::det(dp_dx + I);
    }
};

} // end namespace lido
#endif // __LIDO_SERAC_HELPER_HPP__
#endif // !defined(LIDO_USE_SERAC)

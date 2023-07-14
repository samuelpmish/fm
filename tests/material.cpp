#include "common.hpp"

#include <iomanip>

#include "operations/dot.hpp"
#include "operations/linear_solve.hpp"

#include "operations/random.hpp"

#include "fiber_composite/old_model.hpp"
#include "fiber_composite/new_model.hpp"

#define ANKERL_NANOBENCH_IMPLEMENT
#include "nanobench.h"

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

  ankerl::nanobench::Bench bench;

#if 0
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

  mat3 sigma_original;
  mat3 sigma_new;

  bench.run("original", [&] {
    sigma_original = RotatedFiberComposite3D_Original(
      x, du_dx, alpha0, deltaT, 
      Em, Ez, Ep,
      num, nuzp, nupp,
      Gzp, Vf,
      am, afzz, afpp,
      thermalStress, cylDom
    );
    ankerl::nanobench::doNotOptimizeAway(sigma_original);
  });
#if 0
  auto sigma_new = RotatedFiberComposite3D_New(
    x, du_dx, alpha0, deltaT, 
    Em, Ez, Ep,
    num, nuzp, nupp,
    Gzp, Vf,
    am, afzz, afpp,
    thermalStress, cylDom
  );
#else
  RotatedFiberCompositeModel model({
    Em, Ez, Ep,
    num, nuzp, nupp,
    Gzp, Vf,
    am, afzz, afpp,
    thermalStress, cylDom
  });

  bench.run("new", [&] {
    sigma_new = model(x, du_dx, alpha0, deltaT);
    ankerl::nanobench::doNotOptimizeAway(sigma_new);
  });
#endif

  compare(sigma_original, sigma_new, 1e-13);

}

int main() {
  run_test();
}
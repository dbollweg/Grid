#include <Grid/Grid.h>

using namespace Grid;

template <typename T>
bool has_correct_group_block_structure(const T& U) {
  std::cout << GridLogMessage << "Checking the structure is " << std::endl;
  std::cout << GridLogMessage << "U  =  (   W    X   )  " << std::endl;
  std::cout << GridLogMessage << "      (  -X^*  W^* )  " << std::endl;
  std::cout << GridLogMessage << std::endl;

  const int nsp = Nc / 2;
  Complex i(0., 1.);
  for (int c1 = 0; c1 < nsp; c1++)  // check on W
  {
    for (int c2 = 0; c2 < nsp; c2++) {
      auto W = PeekIndex<ColourIndex>(U, c1, c2);
      auto Wstar = PeekIndex<ColourIndex>(U, c1 + nsp, c2 + nsp);
      auto Ww = conjugate(Wstar);
      auto amizero = sum(W - Ww);
      auto amizeroo = TensorRemove(amizero);
      assert(amizeroo.real() < 10e-6);
      amizeroo *= i;
      assert(amizeroo.real() < 10e-6);
    }
  }

  for (int c1 = 0; c1 < nsp; c1++) {
    for (int c2 = 0; c2 < nsp; c2++) {
      auto X = PeekIndex<ColourIndex>(U, c1, c2 + nsp);
      auto minusXstar = PeekIndex<ColourIndex>(U, c1 + nsp, c2);
      auto minusXx = conjugate(minusXstar);
      auto amizero = sum(X + minusXx);
      auto amizeroo = TensorRemove(amizero);
      assert(amizeroo.real() < 10e-6);
      amizeroo *= i;
      assert(amizeroo.real() < 10e-6);
    }
  }
  return true;
};

template <typename T>
bool is_element_of_sp2n_group(const T& U) {
  LatticeColourMatrixD aux(U.Grid());
  LatticeColourMatrixD identity(U.Grid());
  identity = 1.0;
  LatticeColourMatrixD Omega(U.Grid());
  Sp<Nc>::Omega(Omega);

  std::cout << GridLogMessage << "Check matrix is non-zero " << std::endl;
  assert(norm2(U) > 1e-8);

  std::cout << GridLogMessage << "Unitary check" << std::endl;
  aux = U * adj(U) - identity;
  std::cout << GridLogMessage << "U adjU - 1 = " << norm2(aux) << std::endl;
  assert(norm2(aux) < 1e-8);

  aux = Omega - (U * Omega * transpose(U));
  std::cout << GridLogMessage << "Omega - U Omega transpose(U) = " << norm2(aux)
            << std::endl;
  assert(norm2(aux) < 1e-8);

  std::cout << GridLogMessage
            << "|Det| = " << norm2(Determinant(U)) / U.Grid()->gSites()
            << std::endl;
  assert(norm2(Determinant(U)) / U.Grid()->gSites() - 1 < 1e-8);

  return has_correct_group_block_structure(U);
}

template <typename T>
void test_group_projections(T U) {
  RealD Delta = 666.;
  LatticeColourMatrixD identity(U.Grid());
  identity = 1.0;

  std::cout << GridLogMessage << "#   #   #   #" << std::endl;
  std::cout << GridLogMessage << "Group" << std::endl;
  std::cout << GridLogMessage << "#   #   #   #" << std::endl;
  std::cout << GridLogMessage << std::endl;

  std::string name = "ProjectOnSpGroup";
  std::cout << GridLogMessage << "Testing " << name << std::endl;
  std::cout << GridLogMessage << "Apply to deformed matrix" << std::endl;

  U = U + Delta * identity;
  U = ProjectOnSpGroup(U);
  assert(is_element_of_sp2n_group(U));

  name = "ProjectOnGeneralGroup";
  std::cout << GridLogMessage << "Testing " << name << std::endl;
  std::cout << GridLogMessage << "Apply to deformed matrix" << std::endl;

  U = U + Delta * identity;
  U = Sp<Nc>::ProjectOnGeneralGroup(U);
  assert(is_element_of_sp2n_group(U));

  name = "ProjectOnSpecialGroup";
  std::cout << GridLogMessage << "Testing " << name << std::endl;
  std::cout << GridLogMessage << "Apply to deformed matrix" << std::endl;

  U = U + Delta * identity;
  Sp<Nc>::ProjectOnSpecialGroup(U);
  assert(is_element_of_sp2n_group(U));

  name = "ProjectSpn";
  std::cout << GridLogMessage << "Testing " << name << std::endl;
  std::cout << GridLogMessage << "Apply to deformed matrix" << std::endl;

  U = U + Delta * identity;
  ProjectSpn(U);
  assert(is_element_of_sp2n_group(U));
}

template <typename T>
bool has_correct_algebra_block_structure(const T& U) {
  // this only checks for the anti-hermitian part of the algebra
  const int nsp = Nc / 2;
  Complex i(0., 1.);
  std::cout << GridLogMessage << "Checking the structure is " << std::endl;
  std::cout << GridLogMessage << "U  =  (   W    X   )  " << std::endl;
  std::cout << GridLogMessage << "      (  -X^*  W^* )  " << std::endl;
  for (int c1 = 0; c1 < nsp; c1++)  // check on W
  {
    for (int c2 = 0; c2 < nsp; c2++) {
      auto W = PeekIndex<ColourIndex>(U, c1, c2);
      auto Wstar = PeekIndex<ColourIndex>(U, c1 + nsp, c2 + nsp);
      auto Ww = conjugate(Wstar);
      auto amizero = sum(W - Ww);
      auto amizeroo = TensorRemove(amizero);
      assert(amizeroo.real() < 10e-6);
      amizeroo *= i;
      assert(amizeroo.real() < 10e-6);
    }
  }
  for (int c1 = 0; c1 < nsp; c1++) {
    for (int c2 = 0; c2 < nsp; c2++) {
      auto X = PeekIndex<ColourIndex>(U, c1, c2 + nsp);
      auto minusXstar = PeekIndex<ColourIndex>(U, c1 + nsp, c2);
      auto minusXx = conjugate(minusXstar);
      auto amizero = sum(X + minusXx);
      auto amizeroo = TensorRemove(amizero);
      assert(amizeroo.real() < 10e-6);
      amizeroo *= i;
      assert(amizeroo.real() < 10e-6);
    }
  }

  return true;
}

template <typename T>
bool is_element_of_sp2n_algebra(const T& U) {
  LatticeColourMatrixD aux(U.Grid());
  LatticeColourMatrixD identity(U.Grid());
  identity = 1.0;
  LatticeColourMatrixD Omega(U.Grid());
  Sp<Nc>::Omega(Omega);

  std::cout << GridLogMessage << "Check matrix is non-zero " << std::endl;
  assert(norm2(U) > 1e-8);

  aux = U - adj(U);
  std::cout << GridLogMessage << "T - Tda = " << norm2(aux)
            << " (not supposed to vanish)" << std::endl;

  aux = U + adj(U);
  std::cout << GridLogMessage << "T + Tda = " << norm2(aux)
            << " (supposed to vanish)" << std::endl;
  assert(norm2(aux) - 1 < 1e-8);

  std::cout << GridLogMessage << "Check that Omega T Omega + conj(T) = 0 "
            << std::endl;
  aux = Omega * U * Omega + conjugate(U);
  assert(norm2(aux) < 1e-8);

  return has_correct_algebra_block_structure(U);
}

template <typename T>
void test_algebra_projections(T U) {
  RealD Delta = 666.;
  LatticeColourMatrixD tmp(U.Grid());
  LatticeColourMatrixD identity(U.Grid());
  identity = 1.0;

  std::cout << GridLogMessage << "#   #   #   #" << std::endl;
  std::cout << GridLogMessage << "Algebra" << std::endl;
  std::cout << GridLogMessage << "#   #   #   #" << std::endl;
  std::cout << GridLogMessage << std::endl;

  std::string name = "SpTa";
  std::cout << GridLogMessage << "Testing " << name << std::endl;
  std::cout << GridLogMessage << "Apply to deformed matrix" << std::endl;

  U = U + Delta * identity;
  U = SpTa(U);
  assert(is_element_of_sp2n_algebra(U));

  name = "TaProj";
  std::cout << GridLogMessage << "Testing " << name << std::endl;
  std::cout << GridLogMessage << "Apply to deformed matrix" << std::endl;

  U = U + Delta * identity;
  Sp<Nc>::taProj(U, tmp);
  U = tmp;
  assert(is_element_of_sp2n_algebra(U));
}

int main(int argc, char** argv) {
  Grid_init(&argc, &argv);

  Coordinate latt_size = GridDefaultLatt();
  Coordinate simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  Coordinate mpi_layout = GridDefaultMpi();

  GridCartesian Grid(latt_size, simd_layout, mpi_layout);

  LatticeGaugeField Umu(&Grid);
  LatticeColourMatrixD U(&Grid);

  // Will test resimplectification-related functionalities (from
  // ProjectOnGeneralGroup, ProjectOnSpGroup, ProjectOnSpecialGroup) and projection on the
  // algebra (from SpTa)
  // ProjectOnGeneralGroup, ProjectOnSpGroup project on the non-special group allowi for complex determinants of module 1
  // ProjectOnSpecialGroup projects on the full gauge group providing a determinant equals to 1

  std::vector<int> pseeds({1, 2, 3, 4, 5});
  GridParallelRNG pRNG(&Grid);
  pRNG.SeedFixedIntegers(pseeds);

  SU<Nc>::HotConfiguration(pRNG, Umu);
  U = PeekIndex<LorentzIndex>(Umu, 0);
  test_group_projections(U);
  U = PeekIndex<LorentzIndex>(Umu, 1);
  test_algebra_projections(U);

  Grid_finalize();
}

#include <Grid/Grid.h>

namespace Grid{
  struct ZFParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(ZFParameters,
            int, steps,
            double, step_size,
            int, meas_interval,
            double, maxTau); // for the adaptive algorithm

    
    template <class ReaderClass >
    ZFParameters(Reader<ReaderClass>& Reader){
      read(Reader, "ZeuthenFlow", *this);
    }

  };

  struct ConfParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(ConfParameters,
           std::string, conf_prefix,
            std::string, rng_prefix,
				    int, StartConfiguration,
				    int, EndConfiguration,
            int, Skip);
  
    template <class ReaderClass >
    ConfParameters(Reader<ReaderClass>& Reader){
      read(Reader, "Configurations", *this);
    }

  };


void dfft_solve(GridCartesian& thisgrid, const LatticeFermionD& src, LatticeFermionD &sol, RealD epsilon) {
  FFT theFFT(&thisgrid);
  LatticeFermionD tmp(&thisgrid);
  theFFT.FFT_all_dim(tmp,src,FFT::forward); //tmp = U(omega,0)

  Coordinate latt_size = thisgrid.FullDimensions();

  LatticeComplexD C(&thisgrid);
  LatticeComplexD coor(&thisgrid);

  C=Zero();
  for (int mu = 0; mu < Nd; mu++) {
    RealD PiL = M_PI/latt_size[mu];
    LatticeCoordinate(coor,mu);
    C = C + sin(PiL*coor)*sin(PiL*coor);
  }

  tmp = tmp - 4*epsilon*C*tmp; 
  theFFT.FFT_all_dim(sol,tmp,FFT::backward);
}

}

int main(int argc, char **argv) {
  
  using namespace Grid;

  Grid_init(&argc, &argv);
  GridLogLayout();

  auto latt_size   = GridDefaultLatt();
  auto simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
  auto mpi_layout  = GridDefaultMpi();
  GridCartesian               Grid(latt_size, simd_layout, mpi_layout);
  GridRedBlackCartesian     RBGrid(&Grid);

//  std::vector<int> seeds({1, 2, 3, 4, 5});
  GridSerialRNG sRNG;
  GridParallelRNG pRNG(&Grid);
  //pRNG.SeedFixedIntegers(seeds);

  LatticeGaugeField Umu(&Grid), Uflow(&Grid);
  //SU<Nc>::HotConfiguration(pRNG, Umu);

  typedef Grid::XmlReader       Serialiser;
  Serialiser Reader("input.xml");
  ZFParameters ZFPar(Reader);
  ConfParameters CPar(Reader);
  CheckpointerParameters CPPar(CPar.conf_prefix, CPar.rng_prefix);
  NerscHmcCheckpointer<PeriodicGimplR> CPBin(CPPar);

  for (int conf = CPar.StartConfiguration; conf <= CPar.EndConfiguration; conf+= CPar.Skip){

  CPBin.CheckpointRestore(conf, Umu, sRNG, pRNG);
  SU<Nc>::ColdConfiguration(pRNG,Umu);

  std::cout << std::setprecision(15);
  std::cout << GridLogMessage << "Initial plaquette: "
    << WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu) << std::endl;
  
  //First check: testing on unit gauge, flowed fermion phi should agree with solution from dfft
  LatticeFermionD src(&Grid);
  LatticeFermionD phi(&Grid);
  LatticeFermionD phi_dfft(&Grid);

  random(pRNG, src);

  LatticeFermionD tmp_src(&Grid);
  tmp_src=src;
  for (int i = 0; i < 60; i++) {
    dfft_solve(Grid, tmp_src, phi_dfft, ZFPar.step_size);
    tmp_src=phi_dfft;
  }
  std::cout << GridLogMessage << "Finishd dfft solve!" << std::endl;
  int t=ZFPar.maxTau;
  std::cout << GridLogMessage << "Starting FermionFlow with stepsize " << ZFPar.step_size << std::endl;
  FermionFlow<PeriodicGimplR,ZeuthenAction<PeriodicGimplR>,LatticeFermionD> ZF(ZFPar.step_size, 60, Umu,
					ZFPar.meas_interval,10);

  ZF.smear(Uflow, phi, Umu, src);
  LatticeFermionD diff = phi_dfft - phi;
  
  auto avg = (norm2(PeekIndex<SpinorIndex>(diff,0))+norm2(PeekIndex<SpinorIndex>(diff,1))
             +norm2(PeekIndex<SpinorIndex>(diff,2))+norm2(PeekIndex<SpinorIndex>(diff,3)))/Grid.lSites()/4.0;

  RealD WFlow_plaq = WilsonLoops<PeriodicGimplR>::avgPlaquette(Uflow);
  RealD WFlow_TC   = WilsonLoops<PeriodicGimplR>::TopologicalCharge(Uflow);
  RealD WFlow_T0   = ZF.energyDensityPlaquette(t,Uflow);
  std::cout << GridLogMessage << "Plaquette          "<< conf << "   " << WFlow_plaq << std::endl;
  std::cout << GridLogMessage << "T0                 "<< conf << "   " << WFlow_T0 << std::endl;
  std::cout << GridLogMessage << "TopologicalCharge  "<< conf << "   " << WFlow_TC   << std::endl;

  std::cout << GridLogMessage << "Norm of src = " << norm2(PeekIndex<SpinorIndex>(src,0)) << std::endl;
  std::cout << GridLogMessage << "Norm of phi = " << norm2(PeekIndex<SpinorIndex>(phi,0)) << std::endl;
  std::cout << GridLogMessage << "Norm of phi_dfft = " << norm2(PeekIndex<SpinorIndex>(phi_dfft,0)) << std::endl;
  
  std::cout << GridLogMessage << "Norm of phi_dfft - phi_flow = " << norm2(PeekIndex<SpinorIndex>(diff,0)) << std::endl;
  std::cout << GridLogMessage << "Norm/vol/Ns of phi_dfft - phi_flow = " << avg << std::endl;
  assert(avg < 1e-5);

  //Second check: random gauge transform of unit gauge should agree with gauge transformed solution of dfft
  // and gauge transformed solution of prev check

  LatticeGaugeField Urng(&Grid);
  LatticeGaugeField Uflow_rng(&Grid);
  LatticeFermionD phi_rng(&Grid);

  Urng = Umu;
  LatticeColourMatrix g(&Grid);
  SU<Nc>::RandomGaugeTransform<PeriodicGimplR>(pRNG, Urng, g);

  //rotate src to new gauge;
  src = g * src;
  //rotate phi_dfft to new gauge; 
   
  phi_dfft = g * phi_dfft;

  //rotate phi to new gauge;
  phi = g * phi;

  ZF.smear(Uflow_rng, phi_rng, Urng, src);
  diff = phi_dfft - phi_rng;
  avg = (norm2(PeekIndex<SpinorIndex>(diff,0))+norm2(PeekIndex<SpinorIndex>(diff,1))
        +norm2(PeekIndex<SpinorIndex>(diff,2))+norm2(PeekIndex<SpinorIndex>(diff,3)))/Grid.lSites()/4.0;
  WFlow_plaq = WilsonLoops<PeriodicGimplR>::avgPlaquette(Uflow_rng);
  WFlow_TC   = WilsonLoops<PeriodicGimplR>::TopologicalCharge(Uflow_rng);
  WFlow_T0   = ZF.energyDensityPlaquette(t,Uflow_rng);
  std::cout << GridLogMessage << "Plaquette          "<< conf << "   " << WFlow_plaq << std::endl;
  std::cout << GridLogMessage << "T0                 "<< conf << "   " << WFlow_T0 << std::endl;
  std::cout << GridLogMessage << "TopologicalCharge  "<< conf << "   " << WFlow_TC   << std::endl;

  std::cout << GridLogMessage << "Norm of src = " << norm2(PeekIndex<SpinorIndex>(src,0)) << std::endl;
  std::cout << GridLogMessage << "Norm of phi_rng = " << norm2(PeekIndex<SpinorIndex>(phi_rng,0)) << std::endl;
  std::cout << GridLogMessage << "Norm of phi_dfft = " << norm2(PeekIndex<SpinorIndex>(phi_dfft,0)) << std::endl;
  
  std::cout << GridLogMessage << "Norm of phi_dfft - phi_rng = " << norm2(PeekIndex<SpinorIndex>(diff,0)) << std::endl;
  std::cout << GridLogMessage << "Norm/vol/Ns of phi_dfft - phi_rng = " << avg << std::endl;
  assert(avg < 1e-5);

  diff = phi - phi_rng;
  avg = (norm2(PeekIndex<SpinorIndex>(diff,0))+norm2(PeekIndex<SpinorIndex>(diff,1))
        +norm2(PeekIndex<SpinorIndex>(diff,2))+norm2(PeekIndex<SpinorIndex>(diff,3)))/Grid.lSites()/4.0;
  std::cout << GridLogMessage << "Norm of phi_rng = " << norm2(PeekIndex<SpinorIndex>(phi_rng,0)) << std::endl;
  std::cout << GridLogMessage << "Norm of g*phi = " << norm2(PeekIndex<SpinorIndex>(phi,0)) << std::endl;
  
  std::cout << GridLogMessage << "Norm of g*phi - phi_rng = " << norm2(PeekIndex<SpinorIndex>(diff,0)) << std::endl;
  std::cout << GridLogMessage << "Norm/vol/Ns of phi_dfft - phi_rng = " << avg << std::endl;
  assert(avg < 1e-5);

  LatticeFermionD eta(&Grid);
  random(pRNG, eta);
  LatticeFermionD phi_adj(&Grid);
  auto dotprod_t = innerProduct(eta,phi_rng);
  std::cout << GridLogMessage << "Testing dotprod <eta(t),phi_rng> = " << dotprod_t <<std::endl;
  //Testing adjoint flow
  ZF.smear_adjoint(phi_adj,eta);
  auto dotprod_0 = innerProduct(phi_adj,src);
  std::cout << GridLogMessage << "Norm of eta = " << norm2(PeekIndex<SpinorIndex>(eta,0)) << std::endl;
  std::cout << GridLogMessage << "Testing dotprod <eta(t=0),src> = " << dotprod_0 <<std::endl;
  auto dotdiff = norm(dotprod_t-dotprod_0);
  std::cout << GridLogMessage << "|<eta(t),phi_rng>-<eta(t=0),src>| = " << norm(dotprod_t - dotprod_0) << std::endl;
  assert(dotdiff < 1e-8);
  }
  Grid_finalize();
}  // main

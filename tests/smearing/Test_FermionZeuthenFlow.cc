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


// void fft_solve(GridCartesian& thisgrid, const LatticeFermionD& src, LatticeFermionD &sol, RealD t) {
//   FFT theFFT(&thisgrid);
//   LatticeFermionD tmp(&thisgrid);
//   theFFT.FFT_all_dim(tmp,src,FFT::forward); //tmp = U(omega,0)

//   Coordinate latt_size = thisgrid.FullDimensions();

//   LatticeComplexD C(&thisgrid); //this is going to be omega^2
//   LatticeComplexD coor(&thisgrid);
//   C=Zero();

//   for(int mu =0; mu<4; mu++) {
//     RealD TwoPiL = M_PI * 2.0 / latt_size[mu];
//     LatticeCoordinate(coor,mu);
//     C = C + (TwoPiL*coor)*(TwoPiL*coor); // C = (2pi*n/L)^2 (omega)^2
//   }  
//   tmp = exp(-1.0*C*t) * tmp; //tmp = exp(-omega^2*t)U(omega,0) 

//   theFFT.FFT_all_dim(sol,tmp,FFT::backward);

// }

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
  for (int i = 0; i < 200; i++) {
    dfft_solve(Grid, tmp_src, phi_dfft, ZFPar.step_size);
    tmp_src=phi_dfft;
  }
  std::cout << GridLogMessage << "Finishd dfft solve!" << std::endl;
  int t=ZFPar.maxTau;
  std::cout << GridLogMessage << "Starting FermionFlow with stepsize " << ZFPar.step_size << std::endl;
  FermionFlow<PeriodicGimplR,LatticeFermionD> ZF(ZFPar.step_size, 200,
					ZFPar.meas_interval);

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
  assert(avg < 1e-6);

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
  assert(avg < 1e-6);

  diff = phi - phi_rng;
  avg = (norm2(PeekIndex<SpinorIndex>(diff,0))+norm2(PeekIndex<SpinorIndex>(diff,1))
        +norm2(PeekIndex<SpinorIndex>(diff,2))+norm2(PeekIndex<SpinorIndex>(diff,3)))/Grid.lSites()/4.0;
  std::cout << GridLogMessage << "Norm of phi_rng = " << norm2(PeekIndex<SpinorIndex>(phi_rng,0)) << std::endl;
  std::cout << GridLogMessage << "Norm of g*phi = " << norm2(PeekIndex<SpinorIndex>(phi,0)) << std::endl;
  
  std::cout << GridLogMessage << "Norm of g*phi - phi_rng = " << norm2(PeekIndex<SpinorIndex>(diff,0)) << std::endl;
  std::cout << GridLogMessage << "Norm/vol/Ns of phi_dfft - phi_rng = " << avg << std::endl;
  assert(avg < 1e-6);


  }
  Grid_finalize();
}  // main

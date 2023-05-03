#include <Grid/Grid.h>

namespace Grid{
  struct ZFParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(ZFParameters,
            int, steps,
            double, step_size,
            int, meas_interval,
            double, maxTau); // for the adaptive algorithm

    ZFParameters() {}       

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
}

int main(int argc, char **argv) {
  using namespace Grid;
   ;

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
  Serialiser Reader("input_zf.xml");
  ZFParameters ZFPar(Reader);
  ConfParameters CPar(Reader);
  CheckpointerParameters CPPar(CPar.conf_prefix, CPar.rng_prefix);
  NerscHmcCheckpointer<PeriodicGimplR> CPBin(CPPar);

  for (int conf = CPar.StartConfiguration; conf <= CPar.EndConfiguration; conf+= CPar.Skip){

  CPBin.CheckpointRestore(conf, Umu, sRNG, pRNG);

  std::cout << std::setprecision(15);
  std::cout << GridLogMessage << "Initial plaquette: "
    << WilsonLoops<PeriodicGimplR>::avgPlaquette(Umu) << std::endl;

  int t=ZFPar.maxTau;
  ZeuthenFlow<PeriodicGimplR> ZF(ZFPar.step_size, 20,
					ZFPar.meas_interval);

  ZF.smear(Uflow, Umu);

  RealD WFlow_plaq = WilsonLoops<PeriodicGimplR>::avgPlaquette(Uflow);
  RealD WFlow_TC   = WilsonLoops<PeriodicGimplR>::TopologicalCharge(Uflow);
  RealD WFlow_T0   = ZF.energyDensityPlaquette(t,Uflow);
  std::cout << GridLogMessage << "Plaquette          "<< conf << "   " << WFlow_plaq << std::endl;
  std::cout << GridLogMessage << "T0                 "<< conf << "   " << WFlow_T0 << std::endl;
  std::cout << GridLogMessage << "TopologicalCharge  "<< conf << "   " << WFlow_TC   << std::endl;

  std::cout<< GridLogMessage << " Admissibility check:\n";
  const double sp_adm = 0.067;                // admissible threshold
  const double pl_adm = 1.0-sp_adm/Nc;
  std::cout << GridLogMessage << "   (pl_adm =" << pl_adm << ")\n";

  // Need min and reduce min for this function
  //double sp_max = NC_*(1.0-stpl.plaq_min(U,pl_adm));
  double sp_ave = Nc*(1.0-WFlow_plaq);

  //std::cout<< GridLogMessage << "   sp_max = "        << sp_max <<"\n";
  std::cout<< GridLogMessage << "   sp_ave = "        << sp_ave <<"\n";
  std::cout<< GridLogMessage << "   (sp_admissible = "<< sp_adm <<")\n";
  //std::cout<< GridLogMessage << "   sp_admissible - sp_max = "<<sp_adm-sp_max <<"\n";
  std::cout<< GridLogMessage << "   sp_admissible - sp_ave = "<<sp_adm-sp_ave <<"\n";
  }
  Grid_finalize();
}  // main
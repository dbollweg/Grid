#include <Grid/Grid.h>
#include <Grid/qcd/utils/NprUtils.h>

using namespace std;
using namespace Grid;


int main (int argc, char ** argv)
{
  const int Ls=16; //3GeV
  // const int Ls=12; //4GeV
  Grid_init(&argc,&argv);
  // Double precision grids
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), 
								   GridDefaultSimd(Nd,vComplex::Nsimd()),
								   GridDefaultMpi());
  GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
  GridCartesian         * FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls,UGrid);
  GridRedBlackCartesian * FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls,UGrid);

  std::vector<int> seeds4({1,2,3,4}); 
  GridParallelRNG          RNG4(UGrid);  RNG4.SeedFixedIntegers(seeds4);

  LatticeGaugeField Umu(UGrid);
  LatticeColourMatrix gauge_transformation(UGrid);
  LatticeGaugeField Uprime(UGrid);
 

  std::string config;
  if( argc > 1 && argv[1][0] != '-' )
  {
    std::cout<<GridLogMessage <<"Loading configuration from "<<argv[1]<<std::endl;
    FieldMetaData header;
    NerscIO::readConfiguration(Umu, header, argv[1]);
    config=argv[1];
  }
  else
  {
    std::cout<<GridLogMessage <<"Using hot configuration"<<std::endl;
    //SU<Nc>::HotConfiguration(RNG4,Umu);
    SU<Nc>::ColdConfiguration(Umu);
    config="HotConfig";
  }
  Real alpha=0.06;
  FourierAcceleratedGaugeFixer<PeriodicGimplR>::SteepestDescentGaugeFix(Umu,gauge_transformation,alpha,100000,1.0e-8, 1.0e-8,true,-1);//should be -1
  
  //std::vector<RealD> masses({ 0.000489, 0.0148});//, 0.01483, 0.188});// 4GeV // u/d, s, c 
  std::vector<RealD> masses({0.00046, 0.0224 });//, 0.0186, 0.243}); //3GeV
  int nmass = masses.size();

  std::vector<MobiusFermionR *> FermActs1;
  std::vector<MobiusFermionR *> FermActs2;
  
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"MobiusFermion action as Scaled Shamir kernel"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;

  for(auto mass: masses) {

    RealD M5=1.4;
    RealD b=2;
    RealD c=1;
    MobiusFermionR::ImplParams p1;
    MobiusFermionR::ImplParams p2;
    p1.twist_n_2pi_L = AcceleratorVector<Real,Nd>({-0.41696,0,0.41696,0}); //3GeV twist to p=2GeV
    p2.twist_n_2pi_L = AcceleratorVector<Real,Nd>({0,0.41696,0.41696,0});
    //p1.twist_n_2pi_L = AcceleratorVector<Real,Nd>({0.21277,0,-0.21277,0}); //4GeV twist to p=2GeV
    //p2.twist_n_2pi_L = AcceleratorVector<Real,Nd>({0,-0.21277,-0.21277,0});
 
    // p1.twist_n_2pi_L = AcceleratorVector<Real,Nd>({0.5,0,-0.5,0});
    // p2.twist_n_2pi_L = AcceleratorVector<Real,Nd>({0,-0.5,-0.5,0});
 
    FermActs1.push_back(new MobiusFermionR(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c,p1));
    FermActs2.push_back(new MobiusFermionR(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c,p2));
   
  }
  std::vector<std::pair<Coordinate, Coordinate> > momenta;
 
  for (int i = 2; i < 4; i++) {
  Coordinate mom1({-i,0,i,0});
  Coordinate mom2({0,i,i,0});
   
  std::pair<Coordinate, Coordinate> mompair(mom1,mom2);
  momenta.push_back(mompair);
  }

  NPR<MobiusFermionR> npr_obj_l(*(FermActs1[0]),*(FermActs2[0]));

  std::string nprfilename_ml(config);
  nprfilename_ml += "_npr_res_ml_0.41696twist";
  npr_obj_l.calculate_NPR(momenta, nprfilename_ml);

  NPR<MobiusFermionR> npr_obj_s(*(FermActs1[1]),*(FermActs2[1]));

  std::string nprfilename_ms(config);
  nprfilename_ms += "_npr_res_ms_0.41696twist";
  npr_obj_s.calculate_NPR(momenta, nprfilename_ms);
  
  // NPR<MobiusFermionR> npr_obj_c(*(FermActs1[2]),*(FermActs2[2]));

  // std::string nprfilename_mc(config);
  // nprfilename_mc += "_npr_res_mc_0.5twist";
  // npr_obj_c.calculate_NPR(momenta, nprfilename_mc);

  Grid_finalize();
}
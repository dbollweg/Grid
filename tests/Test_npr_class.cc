#include <Grid/Grid.h>
#include <Grid/qcd/utils/NprUtils.h>

using namespace std;
using namespace Grid;


int main (int argc, char ** argv)
{
  const int Ls=12;

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
  
  std::vector<RealD> masses({ 0.00035, 0.0142, 0.188});//,0.04,0.45} ); // u/d, s, c 

  int nmass = masses.size();

  std::vector<MobiusFermionR *> FermActs;
  
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"MobiusFermion action as Scaled Shamir kernel"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;

  for(auto mass: masses) {

    RealD M5=1.4;
    RealD b=2;
    RealD c=1;
    
    FermActs.push_back(new MobiusFermionR(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c));
   
  }
  std::vector<std::pair<Coordinate, Coordinate> > momenta;
 
  for (int i = 2; i < 4; i++) {
  Coordinate mom1({-i,0,i,0});
  Coordinate mom2({0,i,i,0});
   
  std::pair<Coordinate, Coordinate> mompair(mom1,mom2);
  momenta.push_back(mompair);
  }

  NPR<MobiusFermionR> npr_obj_l(*(FermActs[0]));

  std::string nprfilename_ml(config);
  nprfilename_ml += "_npr_res_ml";
  npr_obj_l.calculate_NPR(momenta, nprfilename_ml);

  NPR<MobiusFermionR> npr_obj_s(*(FermActs[1]));

  std::string nprfilename_ms(config);
  nprfilename_ms += "_npr_res_ms";
  npr_obj_s.calculate_NPR(momenta, nprfilename_ms);
  
  NPR<MobiusFermionR> npr_obj_c(*(FermActs[2]));

  std::string nprfilename_mc(config);
  nprfilename_mc += "_npr_res_mc";
  npr_obj_c.calculate_NPR(momenta, nprfilename_mc);

  Grid_finalize();
}
#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

void MakePhase(Coordinate &mom,LatticeComplex &phase)
{
  GridBase *grid = phase.Grid();
  auto latt_size = grid->GlobalDimensions();
  ComplexD ci(0.0,1.0);
  phase=Zero();

  LatticeComplex coor(phase.Grid());
  for(int mu=0;mu<Nd;mu++){
    RealD TwoPiL =  M_PI * 2.0/ latt_size[mu];
    LatticeCoordinate(coor,mu);
    phase = phase + (TwoPiL * mom[mu]) * coor;
  }
  phase = exp(phase*ci);
}


void VolumeSource(Coordinate &mom, LatticePropagator &source)
{
  source=1.0;
  
  LatticeComplex phase(source.Grid());
  MakePhase(mom,phase);

  source = phase * source;
}


template<class Action>
void Solve(Action &D,LatticePropagator &source,LatticePropagator &propagator)
{
  GridBase *UGrid = D.GaugeGrid();
  GridBase *FGrid = D.FermionGrid();

  LatticeFermion src4  (UGrid); 
  LatticeFermion src5  (FGrid); 
  LatticeFermion result5(FGrid);
  LatticeFermion result4(UGrid);
  
  ConjugateGradient<LatticeFermion> CG(1.0e-8,100000);
  SchurRedBlackDiagMooeeSolve<LatticeFermion> schur(CG);
  ZeroGuesser<LatticeFermion> ZG; 
  for(int s=0;s<Nd;s++){
    for(int c=0;c<Nc;c++){
      PropToFerm<Action>(src4,source,s,c);

      D.ImportPhysicalFermionSource(src4,src5);

      result5=Zero();
      schur(D,src5,result5,ZG);
      std::cout<<GridLogMessage
	       <<"spin "<<s<<" color "<<c
	       <<" norm2(src5d) "   <<norm2(src5)
               <<" norm2(result5d) "<<norm2(result5)<<std::endl;

      D.ExportPhysicalFermionSolution(result5,result4);

      FermToProp<Action>(propagator,result4,s,c);
    }
  }
}

template<class Action>
LatticePropagator PhasedPropagator(Action &D, Coordinate p_in, LatticeColourMatrix &gauge_transformation)
{
  //Create Volumesource, rotate to landau gauge, solve for propagator, multiply phase
  GridBase *UGrid = D.GaugeGrid();

  LatticePropagator src4 (UGrid);
  LatticePropagator result (UGrid);
  LatticeComplex phase (UGrid);

  VolumeSource(p_in, src4);

  src4 = gauge_transformation * src4;

  Solve(D, src4, result);
  
  Coordinate minus_p_in(p_in);
  for(size_t i =0;i < p_in.size() ;i++){
    minus_p_in[i] = -1.0 * p_in[i];
  }
  MakePhase(minus_p_in, phase); //Multiply with exp(-i p_in x) (see equation (8) in https://arxiv.org/abs/1006.0422v2)
  result = phase * result;
  return result;
}

template<class Action>
auto ExternalLeg(Action &D, Coordinate p, LatticeColourMatrix &gauge_transformation)
{
  GridBase *UGrid = D.GaugeGrid();

  LatticePropagator G(UGrid);

  G = PhasedPropagator(D, p, gauge_transformation);
  return sum(G);
}

template<class Action>
auto BilinearVertex(Action &D, Coordinate p1, Coordinate p2, LatticeColourMatrix &gauge_transformation)
{
  //Create Phased propagators with momenta p1 and p2, then construct vertices for the 16 Gamma combinations
  GridBase *UGrid = D.GaugeGrid();
  
  LatticePropagator G1(UGrid);
  LatticePropagator G2(UGrid);

  G1 = PhasedPropagator(D, p1, gauge_transformation);
  G2 = PhasedPropagator(D, p2, gauge_transformation);

  std::array<SpinColourMatrix, 16> vertex;

  
  for (size_t i = 0; i < 16; i++) 
  {
    Gamma G5(Gamma::Algebra::Gamma5);
    Gamma Gi(Gamma::gall[i]);
    vertex[i] = sum((G5 * adj(G1) * G5) * Gi * G2);
  }

  return vertex;

}

template<class Action>
auto FourQuarkOperator(Action &D, Coordinate p1, Coordinate p2, LatticeColourMatrix &gauge_transformation)
{
  GridBase *UGrid = D.GaugeGrid();

  LatticePropagator G1(UGrid);
  LatticePropagator G2(UGrid);

  G1 = PhasedPropagator(D, p1, gauge_transformation);
  G2 = PhasedPropagator(D, p2, gauge_transformation);

  std::array<SpinColourMatrix, 16> FourQuarkOp;

  for (size_t i = 0; i < 16; i++) 
  {
    Gamma G5(Gamma::Algebra::Gamma5);
    Gamma Gi(Gamma::gall[i]);
    FourQuarkOp[i] = sum(((G5 * adj(G1) * G5) * Gi * G2) * ((G5 * adj(G1) * G5) * Gi * G2));
  }

  return FourQuarkOp;
}

int main (int argc, char ** argv)
{
  const int Ls=8;

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

  Real alpha=0.1;

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
    SU<Nc>::HotConfiguration(RNG4,Umu);
    config="HotConfig";
  }

  std::vector<RealD> masses({ 0.03});//,0.04,0.45} ); // u/d, s, c 

  int nmass = masses.size();

  std::vector<MobiusFermionR *> FermActs;
  
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"MobiusFermion action as Scaled Shamir kernel"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;

  for(auto mass: masses) {

    RealD M5=1.0;
    RealD b=1.5;// Scale factor b+c=2, b-c=1
    RealD c=0.5;
    
    FermActs.push_back(new MobiusFermionR(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c));
   
  }
  
  FourierAcceleratedGaugeFixer<PeriodicGimplR>::SteepestDescentGaugeFix(Umu,gauge_transformation,alpha,10000,1.0e-12, 1.0e-12,true,Nd);
  Uprime = Umu;
  SU<Nc>::GaugeTransform(Uprime, gauge_transformation);
  

  //  Coordinate mom1({-1,0,1,0});
  //   Coordinate mom2({0,1,1,0});

  std::vector<std::array<SpinColourMatrix, 16> > bilinear_vertices; 
  std::vector<SpinColourMatrix> leg_1; 
  std::vector<SpinColourMatrix> leg_2; 

  for(int m=0;m<nmass;m++) {
    for (int mom_factor = 1; mom_factor < 5; mom_factor++) 
    {
      Coordinate mom1({-mom_factor,0,mom_factor,0});
      Coordinate mom2({0,mom_factor,mom_factor,0});
      std::cout << GridLogMessage << "momentum factor " << mom_factor << std::endl;
      std::cout << GridLogMessage << "Calculating ExternalLeg 1" << std::endl;
      leg_1.push_back(ExternalLeg(*(FermActs[m]), mom1, gauge_transformation));

      std::cout << GridLogMessage << "Calculating ExternalLeg 2" << std::endl;
      leg_2.push_back(ExternalLeg(*(FermActs[m]), mom2, gauge_transformation));

      std::cout << GridLogMessage << "Calculating Bilinear vertices" << std::endl;
      bilinear_vertices.push_back(BilinearVertex(*(FermActs[m]),mom1, mom2, gauge_transformation));
    }
  }

  for(int m=0;m<nmass;m++) {
    for (int mom_factor = 1; mom_factor < 5; mom_factor++) 
    {
      int index = m*4+mom_factor-1;
      std::cout << GridLogMessage << "m = " << masses[m] << ", mom_factor = " << mom_factor << std::endl;
      std::cout << GridLogMessage << "leg 1:" << std::endl << leg_1[index] <<std::endl << std::endl;
      std::cout << GridLogMessage << "leg 2:" << std::endl << leg_2[index] << std::endl << std::endl;
      std::cout << GridLogMessage << "unamputated vertex [gamma5]:" << std::endl << bilinear_vertices[index][1] << std::endl << std::endl;
     
    }
  }

  
  Grid_finalize();
}




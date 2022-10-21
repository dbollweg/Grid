#include <Grid/Grid.h>

using namespace std;
using namespace Grid;


SpinColourMatrix invert_eigen(SpinColourMatrix input)
{
    SpinColourMatrix output;
    Eigen::MatrixXcd mat(12,12);
    for (int i1 = 0; i1 < Nd; i1++) {
        for (int j1 = 0; j1 < Nc; j1++) {
            for (int i2 = 0; i2 < Nd; i2++) {
                for (int j2 = 0; j2 < Nc; j2++) {
                    int index_1 = i1*Nc+j1;
                    int index_2 = i2*Nc+j2;

                    mat(index_1,index_2) = input()(i1,i2)(j1,j2);
                }
            }
        }
    }

    Eigen::MatrixXcd invmat(12,12);
    invmat = mat.inverse();

    for (int i1 = 0; i1 < Nd; i1++) {
        for (int j1 = 0; j1 < Nc; j1++) {
            for (int i2 = 0; i2 < Nd; i2++) {
                for (int j2 = 0; j2 < Nc; j2++) {
                    int index_1 = i1*Nc+j1;
                    int index_2 = i2*Nc+j2;

                    output()(i1,i2)(j1,j2) = invmat(index_1,index_2);
                }
            }
        }
    }
    return output;
}

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
 // std::cout<< GridLogMessage << "Testing VolumeSource" << source <<std::endl;
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
  ZeroGuesser<LatticeFermion> ZG; // Could be a DeflatedGuesser if have eigenvectors
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

      std::cout << GridLogMessage << "spin " << s << " color " << c << " norm2(src4) " << norm2(src4) << " norm2(result4) " << norm2(result4) <<std::endl;
      FermToProp<Action>(propagator,result4,s,c);
    }
  }
}

template<class Action>
void Dslash_mult(Action &D, LatticePropagator &source, LatticePropagator &result)
{
  GridBase *UGrid = D.GaugeGrid();
  GridBase *FGrid = D.FermionGrid();

  LatticeFermion src4 (UGrid);
  LatticeFermion src5 (FGrid);
  LatticeFermion result5(FGrid);
  LatticeFermion result4(UGrid);

  for(int s=0;s<Nd;s++){
    for(int c=0;c<Nc;c++) {
      PropToFerm<Action>(src4,source,s,c);
      D.ImportPhysicalFermionSource(src4,src5);

      std::cout << GridLogMessage << "spin " << s << " color " << c << " norm2(src5) " << norm2(src5) << " norm2(result5) " << norm2(result5) <<std::endl;
      D.M(src5,result5);
      std::cout << GridLogMessage << "spin " << s << " color " << c << " norm2(src5) " << norm2(src5) << " norm2(result5) " << norm2(result5) <<std::endl;
      D.ExportPhysicalFermionSource(result5,result4);

      std::cout << GridLogMessage << "spin " << s << " color " << c << " norm2(src4) " << norm2(src4) << " norm2(result4) " << norm2(result4) <<std::endl;
      FermToProp<Action>(result,result4,s,c);
    }
  }
}

class NprFile: Serializable {
  public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(NprFile, std::vector<SpinColourMatrix>, data);
};

template<class action>
class NPR {
    public:
    NPR(action &D, LatticeGaugeField &Umu) : _D(D), _UGrid(_D.GaugeGrid()), _G1(_UGrid), _G2(_UGrid), _gauge_transformation(_UGrid) 
    {
        FourierAcceleratedGaugeFixer<PeriodicGimplR>::SteepestDescentGaugeFix(Umu,_gauge_transformation,_alpha,10000,1.0e-12, 1.0e-12,true,-1);//should be -1    
    }
    
    LatticePropagator PhasedPropagator(Coordinate p);

    auto BilinearVertex(LatticePropagator G1, LatticePropagator G2) 
    {
        std::array<SpinColourMatrix, 16> vertex;
        
    
        for (size_t i = 0; i < 16; i++) 
        {
            Gamma G5(Gamma::Algebra::Gamma5);
            Gamma Gi(Gamma::gall[i]);
            vertex[i] = sum((G5 * adj(G1) * G5) * Gi * adj(G2));
        }

        return vertex;
    };

    SpinColourMatrix ExternalLeg(LatticePropagator G) {return sum(G);};

    void calculate_NPR(std::vector<std::pair<Coordinate, Coordinate> > momenta);

    void save_results();

    private:
    action &_D;
    GridBase *_UGrid;
    LatticePropagator _G1;
    LatticePropagator _G2;
    LatticeColourMatrix _gauge_transformation;
    RealD _alpha = 0.1;

};

template<class action>
LatticePropagator NPR<action>::PhasedPropagator(Coordinate p)
{
    LatticePropagator src(_UGrid);
    LatticePropagator result(_UGrid);
    LatticePropagator tester(_UGrid);
    
    LatticeComplex phase(_UGrid);
    MakePhase(p,phase);
    VolumeSource(p, src);
    Solve(_D, src, result); 

    result = adj(phase) * result;

    return result;
};

template <class action>
void NPR<action>::calculate_NPR(std::vector<std::pair<Coordinate, Coordinate> > momenta) 
{
    for (auto mom: momenta) {

        _G1 = PhasedPropagator(mom.first);
        _G2 = PhasedPropagator(mom.second);
        NprFile NF;
        
        std::array<SpinColourMatrix, 16> bilinear_vertices = BilinearVertex(_G1,_G2);
        SpinColourMatrix leg1 = ExternalLeg(_G1);
        SpinColourMatrix leg2 = ExternalLeg(_G2);

        SpinColourMatrix leginv1 = invert_eigen(leg1);
        SpinColourMatrix leginv2 = invert_eigen(leg2);
        LatticePropagator test(_UGrid);
        test = 1.0;

        //This part is for testing on unit config, should give 12 for each gamma

        std::cout << GridLogMessage << "test inversion 1 " << leg1*leginv1 << std::endl << std::endl;
        std::cout << GridLogMessage << "printing leg1 " << leg1 << std::endl << std::endl;
        std::cout << GridLogMessage << "printing invleg1 " << leginv1 << std::endl << std::endl;
        
        std::cout << GridLogMessage << "trace invleg1*leg1 " << TensorRemove(trace(leg1*leginv1)) << std::endl;
        std::cout << GridLogMessage << "testing norm of sum" << TensorRemove(trace(sum(test))) << std::endl;

        std::cout << GridLogMessage << "test inversion 2" << leg2*leginv2 << std::endl << std::endl;
        std::cout << GridLogMessage << "printing leg2 " << leg2 << std::endl << std::endl;
        std::cout << GridLogMessage << "printing invleg2 " << leginv2 << std::endl << std::endl;
        
        std::cout << GridLogMessage << "trace invleg2*leg2 " << TensorRemove(trace(leg2*leginv2)) << std::endl;
        for (size_t i = 0; i < 16; i++) 
        {

            Gamma G(Gamma::gall[i]);
            Gamma G5(Gamma::Algebra::Gamma5);
            auto tmp = leginv1 * bilinear_vertices[i] * (G5 * adj(leginv2) * G5);
            NF.data.push_back(tmp);
            auto vol=_UGrid->_gsites;
            RealD volD = vol;
            auto result = TensorRemove(trace(tmp * G));
            std::cout << GridLogMessage << "test gamma " << i << " " << result << std::endl;
            
        }


        //testing I/O 
        std::string file = "test_mom_" + std::to_string(mom.first[0]) + std::to_string(mom.first[1]) + std::to_string(mom.first[2]) + std::to_string(mom.first[3]) + ".txt";
        
        TextWriter WR(file); 
        write(WR,"NprFile",NF);
        
        TextReader RD(file);
        read(RD, "NprFile",NF);
        std::vector<SpinColourMatrix> readtest;
        readtest = NF.data;

        for (auto dat: readtest) {
          std::cout << GridLogMessage << "Testing FileReader: " << dat << std::endl;
        }

        for (auto dat: NF.data) {
          std::cout << GridLogMessage << "Testing FileRead (reference data) " << dat << std::endl;
        }

    }
};

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
    //SU<Nc>::HotConfiguration(RNG4,Umu);
    SU<Nc>::ColdConfiguration(Umu);
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
  std::vector<std::pair<Coordinate, Coordinate> > momenta;
  Coordinate mom1({-1,0,1,0});
  Coordinate mom2({0,1,1,0});
  Coordinate mom0({0,0,0,0});
  std::pair<Coordinate, Coordinate> pair1(mom0,mom0);
  momenta.push_back(pair1);

  NPR<MobiusFermionR> npr_obj(*(FermActs[0]), Umu);
  npr_obj.calculate_NPR(momenta);





  
  Grid_finalize();
}
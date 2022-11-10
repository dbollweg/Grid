#include <Grid/Grid.h>

using namespace std;
using namespace Grid;

ComplexD mean(const std::vector<ComplexD>& data)
{
    int N = data.size();
    ComplexD mean(0.0);
    for(int i=0; i<N; ++i){ mean += data[i]; }
    return mean/ComplexD(N);
}

ComplexD jack_mean(const std::vector<ComplexD>& data, int sample)
{
    int N = data.size();
    ComplexD mean(0.0);
    for(int i=0; i<N; ++i){ if(i != sample){ mean += data[i]; } }
    return mean/ComplexD(N-1);
}

RealD jack_std(const std::vector<ComplexD>& jacks, ComplexD mean)
{
    int N = jacks.size();
    ComplexD std(0.0);
    for(int i=0; i<N; ++i){ std += (jacks[i]-mean)*adj(jacks[i]-mean); }
    return std::sqrt(RealD(N-1)/RealD(N)*std.real());
}

std::vector<ComplexD> jack_stats(const std::vector<ComplexD>& data)
{
    int N = data.size();
    std::vector<ComplexD> jack_samples(N);
    std::vector<ComplexD> jack_stats(2);

    jack_stats[0] = mean(data);
    for(int i=0; i<N; i++){ jack_samples[i] = jack_mean(data,i); }
    jack_stats[1] = ComplexD(jack_std(jack_samples, jack_stats[0]),0);
    return jack_stats;
}

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


class NprFile: Serializable {
  public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(NprFile, std::vector<SpinColourMatrix>, data_bv, SpinColourMatrix, data_leg1, SpinColourMatrix, data_leg2);
};



//TODO put together stats analysis + contraction class
template<class action>
class NPR_analyze {
  public:
  NPR_analyze(action &D, std::vector<std::string> filenames) : _D(D), _UGrid(_D.GaugeGrid()), _filenames(filenames) {
    
    for (size_t i=0; i < _filenames.size(); i++) {

      NprFile tmpfile;
      data.push_back(tmpfile);
      BinaryReader BR(_filenames[i]);
      
      try{ read(BR,"NprFile",data[i]); }
      catch(const std::bad_alloc&) {
        std::cout << GridLogError << "Error loading " << _filenames[i] << std::endl;
      }
    
    }
  
  }

  auto Average();


  private:
  action &_D;
  GridBase *_UGrid;
  std::vector<std::string> _filenames;
  std::vector<NprFile> data;
};

template<class action>
auto NPR_analyze<action>::Average() {

    std::vector<ComplexD> tmp_data_leg1(_filenames.size()); 
    std::vector<ComplexD> tmp_data_leg2(_filenames.size());
    std::array<std::vector<ComplexD>,16> tmp_data_bv;
    for (size_t g = 0; g < 16; g++) {
      tmp_data_bv[g].resize(_filenames.size());
    }

    SpinColourMatrix avg_leg1;
    SpinColourMatrix avg_leg2;
    std::array<SpinColourMatrix,16> avg_bv;   


    SpinColourMatrix err_leg1;
    SpinColourMatrix err_leg2;
    std::array<SpinColourMatrix,16> err_bv; 

      for (int s1 = 0; s1 < Nd; s1++) {
        for (int c1 = 0; c1 < Nc; c1++) {
          for (int s2 = 0; s2 < Nd; s2++) {
            for (int c2 = 0; c2 < Nc; c2++) {


              for (size_t i = 0; i < data.size(); i++) {  

                tmp_data_leg1[i]=data[i].data_leg1()(s1,s2)(c1,c2);
                tmp_data_leg2[i]=data[i].data_leg2()(s1,s2)(c1,c2);
                for (size_t g = 0; g < 16; g++) { 
                  tmp_data_bv[g][i]=data[i].data_bv[g]()(s1,s2)(c1,c1);
                }

            }


            std::vector<ComplexD> leg1 = jack_stats(tmp_data_leg1);
            std::vector<ComplexD> leg2 = jack_stats(tmp_data_leg2);
            std::array<std::vector<ComplexD>,16> bv;
            for (size_t g = 0; g < 16; g++) {
              bv[g] = jack_stats(tmp_data_bv[g]);
            }

            avg_leg1()(s1,s2)(c1,c2) = leg1[0];
            avg_leg2()(s1,s2)(c1,c2) = leg2[0];

            for (size_t g = 0; g < 16; g++) {
              avg_bv[g]()(s1,s2)(c1,c2) = bv[g][0];
            }


            err_leg1()(s1,s2)(c1,c2) = leg1[1];
            err_leg2()(s1,s2)(c1,c2) = leg2[1];

            for (size_t g = 0; g < 16; g++) {
              err_bv[g]()(s1,s2)(c1,c2) = bv[g][1];
            }
            
          }          
        }
      }
    }

  //do contraction now (How to do error propagation?)


  SpinColourMatrix leginv1 = invert_eigen(avg_leg1);
  SpinColourMatrix leginv2 = invert_eigen(avg_leg2);
  std::array<ComplexD, 16> result;
  for (size_t i = 0; i < 16; i++) 
        {

            Gamma G(Gamma::gall[i]);
            Gamma G5(Gamma::Algebra::Gamma5);
            auto tmp = (G5 * adj(leginv1) * G5)* avg_bv[i] * leginv2;
           
            auto vol=_UGrid->_gsites;
            RealD volD = vol;
            
            result[i] = volD*TensorRemove(trace(tmp * G))/12.0;
            std::cout << GridLogMessage << "Lambda(gamma) " << i << " " << result[i] << std::endl;
            
        }
  return result;
}

template<class action>
class NPR {
    public:
    NPR(action &D) : _D(D), _UGrid(_D.GaugeGrid()), _G1(_UGrid), _G2(_UGrid)
    {
            
        
    }
    
    LatticePropagator PhasedPropagator(Coordinate p);

    auto BilinearVertex(LatticePropagator G1, LatticePropagator G2) 
    {
        std::array<SpinColourMatrix, 16> vertex;
        
    
        for (size_t i = 0; i < 16; i++) 
        {
            Gamma G5(Gamma::Algebra::Gamma5);
            Gamma Gi(Gamma::gall[i]);
            vertex[i] = sum((G5 * adj(G1) * G5) * Gi * G2);
        }

        return vertex;
    };

    SpinColourMatrix ExternalLeg(LatticePropagator G) {return sum(G);};

    //for testing only
    void test_NPR(std::vector<std::pair<Coordinate, Coordinate> > momenta);


    void calculate_NPR(std::vector<std::pair<Coordinate, Coordinate> > momenta, std::string base_filename);


    private:
    action &_D;
    GridBase *_UGrid;
    LatticePropagator _G1;
    LatticePropagator _G2;
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
void NPR<action>::test_NPR(std::vector<std::pair<Coordinate, Coordinate> > momenta) 
{
    for (auto mom: momenta) {

        _G1 = PhasedPropagator(mom.first);
        _G2 = PhasedPropagator(mom.second);
        
        
        std::array<SpinColourMatrix, 16> bilinear_vertices = BilinearVertex(_G1,_G2);
        SpinColourMatrix leg1 = ExternalLeg(_G1);
        SpinColourMatrix leg2 = ExternalLeg(_G2);

        SpinColourMatrix leginv1 = invert_eigen(leg1);
        SpinColourMatrix leginv2 = invert_eigen(leg2);
        LatticePropagator test(_UGrid);

        for (size_t i = 0; i < 16; i++) 
        {

            Gamma G(Gamma::gall[i]);
            Gamma G5(Gamma::Algebra::Gamma5);
            auto tmp = ((G5 * adj(leginv1) * G5) * bilinear_vertices[i] * leginv2);
            std::cout << GridLogMessage << "Pi " << i << " " << tmp << std::endl;
            auto vol=_UGrid->_gsites;
            RealD volD = vol;
            auto result = volD*TensorRemove(trace(tmp * G))/12.0;
            std::cout << GridLogMessage << "test gamma " << i << " " << result << std::endl;
            
        }


        

    }
};




template <class action>
void NPR<action>::calculate_NPR(std::vector<std::pair<Coordinate, Coordinate> > momenta, std::string base_filename) 
{
    for (auto mom: momenta) {
        NprFile save_file;

        std::string filename = base_filename + "_p1_";
        for (auto m1: mom.first){
          filename += std::to_string(m1);
        }

        filename += "_p2_";
        for (auto m2: mom.second){
          filename += std::to_string(m2);
        }

        filename += ".dat";
        
        std::cout << GridLogMessage << "first phased propagator" << std::endl;
        _G1 = PhasedPropagator(mom.first);
        std::cout << GridLogMessage << "second phased propagator" << std::endl;
        _G2 = PhasedPropagator(mom.second);

        std::cout << GridLogMessage << "Instantiating BilinearVertex" << std::endl;
        std::array<SpinColourMatrix, 16> bilinear_vertices = BilinearVertex(_G1,_G2);

        std::cout << GridLogMessage << "First external leg" <<std::endl;
        SpinColourMatrix leg1 = ExternalLeg(_G1);

        std::cout << GridLogMessage << "Second external leg" <<std::endl;
        SpinColourMatrix leg2 = ExternalLeg(_G2);

        SpinColourMatrix leginv1 = invert_eigen(leg1);
        SpinColourMatrix leginv2 = invert_eigen(leg2);

        for (size_t i = 0; i < 16; i++) 
        {

            Gamma G(Gamma::gall[i]);
            Gamma G5(Gamma::Algebra::Gamma5);
            auto tmp = (G5 * adj(leginv1) * G5) * bilinear_vertices[i] * leginv2;
           
            auto vol=_UGrid->_gsites;
            RealD volD = vol;
            auto result = volD*TensorRemove(trace(tmp * G))/12.0;
            std::cout << GridLogMessage << "test gamma " << i << " " << result << std::endl;
            
        }


        BinaryWriter BWR(filename);
        for (auto bv: bilinear_vertices) {
          save_file.data_bv.push_back(bv);
        }
        save_file.data_leg1 = leg1;
        save_file.data_leg2 = leg2;

        write(BWR,"NprFile", save_file);

    }
};


int main (int argc, char ** argv)
{
  const int Ls=16;

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
  Real alpha=0.1;
  FourierAcceleratedGaugeFixer<PeriodicGimplR>::SteepestDescentGaugeFix(Umu,gauge_transformation,alpha,100000,1.0e-10, 1.0e-10,false,-1);//should be -1
  
  std::vector<RealD> masses({ 0.01});//,0.04,0.45} ); // u/d, s, c 

  int nmass = masses.size();

  std::vector<MobiusFermionR *> FermActs;
  
  std::cout<<GridLogMessage <<"======================"<<std::endl;
  std::cout<<GridLogMessage <<"MobiusFermion action as Scaled Shamir kernel"<<std::endl;
  std::cout<<GridLogMessage <<"======================"<<std::endl;

  for(auto mass: masses) {

    RealD M5=1.8;
    RealD b=1.5;
    RealD c=0.5;
    
    FermActs.push_back(new MobiusFermionR(Umu,*FGrid,*FrbGrid,*UGrid,*UrbGrid,mass,M5,b,c));
   
  }
  std::vector<std::pair<Coordinate, Coordinate> > momenta;
 
  //Coordinate mom0({0,0,0,0});
    Coordinate mom1({-2,0,2,0});
    Coordinate mom2({0,2,2,0});
   
    std::pair<Coordinate, Coordinate> mompair(mom1,mom2);
    momenta.push_back(mompair);
  

  NPR<MobiusFermionR> npr_obj(*(FermActs[0]));
  //npr_obj.test_NPR(momenta);

  std::string nprfilename(config);
  nprfilename += "_npr_res";
  npr_obj.calculate_NPR(momenta, nprfilename);

  /*
  testfilename += "_p1_0000_p2_0000.dat";
  std::vector<std::string> filenames;
  filenames.push_back(testfilename);
  filenames.push_back("testfile1_p1_0000_p2_0000.dat");
  filenames.push_back("testfile2_p1_0000_p2_0000.dat");

  NPR_analyze<MobiusFermionR> npr_reader(*(FermActs[0]),filenames);

  auto tmp = npr_reader.Average();*/
  
  Grid_finalize();
}
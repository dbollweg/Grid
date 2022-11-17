#pragma once

NAMESPACE_BEGIN(Grid);

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
class NPR_analyze {
  public:
  NPR_analyze(std::vector<std::string> filenames, Coordinate lattDim) : _filenames(filenames), vol(lattDim[0]*lattDim[1]*lattDim[2]*lattDim[3]) {
    
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

  void Average();
  template<class T> 
  std::vector<T> CreateBlocks(std::vector<T> &input_data);


  void Test();
  private:
  std::vector<std::string> _filenames;
  std::vector<NprFile> data;
  int vol;
};

void NPR_analyze::Test() {
    int N = data.size();
    std::vector<SpinColourMatrix> data_leg1(N);
    std::vector<SpinColourMatrix> jack_means_leg1(N);

    for (int i = 0; i < N; i++) {
        data_leg1[i]=data[i].data_leg1;
         std::cout << GridLogMessage << "Printing data of leg1 index " << i << "   " << data_leg1[i] << std::endl;
    }

    jack_means_leg1=CreateBlocks(data_leg1);


    for (int i = 0; i < N; i++) {
        std::cout << GridLogMessage << "Testing jack_means on leg1 index " << i << "   " << jack_means_leg1[i] << std::endl;
    }


}

template<class T>
std::vector<T> NPR_analyze::CreateBlocks(std::vector<T> &input_data) {
    int N = input_data.size();

    std::vector<T> blocks_data(N);     


    for (int sample = 0; sample < N; ++sample) {
        T mean(0.0);

        for(int i=0; i<N; ++i){ if(i != sample){ blocks_data[sample] += input_data[i]; } }

        blocks_data[sample] = (1.0/(N-1)) * blocks_data[sample];
    }

    return blocks_data;
}

void NPR_analyze::Average() {
  int N = data.size();

  std::vector<SpinColourMatrix> leg1_blocks(N);
  std::vector<SpinColourMatrix> leg2_blocks(N);
  std::array<std::vector<SpinColourMatrix>, 16> bv_blocks;
  for (size_t i = 0; i < 16; i++) {
          
    bv_blocks[i].resize(N);
  } 
  for (int j = 0; j < N; j++) {
    leg1_blocks[j]=data[j].data_leg1;
    leg2_blocks[j]=data[j].data_leg2;
    
    for (size_t i = 0; i < 16; i++) 
    {
      bv_blocks[i][j]=data[j].data_bv[i];
    }
  }

  leg1_blocks = CreateBlocks(leg1_blocks);
  leg2_blocks = CreateBlocks(leg2_blocks);

  for (size_t i = 0; i < 16; i++) {
    bv_blocks[i] = CreateBlocks(bv_blocks[i]);
  }

  std::vector<SpinColourMatrix> inv_leg1(N);
  std::vector<SpinColourMatrix> inv_leg2(N);
  
  for (int j = 0; j < N; j++) {
    inv_leg1[j]=invert_eigen(leg1_blocks[j]);
    inv_leg2[j]=invert_eigen(leg2_blocks[j]);
  }

  
  
  
  std::array<std::vector<ComplexD>,16> block_results;

  for (size_t i = 0; i < 16; i++) {
          
    block_results[i].resize(N);
  } 

  for (int j = 0; j < N; j++) {
  for (size_t i = 0; i < 16; i++) 
        {

            Gamma G(Gamma::gall[i]);
            Gamma G5(Gamma::Algebra::Gamma5);
            auto tmp = (G5 * adj(inv_leg1[j]) * G5)* bv_blocks[i][j] * inv_leg2[j];

            block_results[i][j] = RealD(vol)*TensorRemove(trace(tmp * G))/12.0;
         
        }
  }

  for (size_t i = 0; i < 16; i++) {
    std::vector<ComplexD> result(2);
    result = jack_stats(block_results[i]);
    std::cout << GridLogMessage << "Lambda(gamma) " << i << ": " << result[0] << " +/- " << result[1]<<std::endl;
  }
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
}

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
}




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

        std::array<SpinColourMatrix, 16> bilinear_vertices = BilinearVertex(_G1,_G2);

        SpinColourMatrix leg1 = ExternalLeg(_G1);

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
            std::cout << GridLogMessage << "Vertex (gamma_" << i << " )=" << result << std::endl;
            
        }


        BinaryWriter BWR(filename);
        for (auto bv: bilinear_vertices) {
          save_file.data_bv.push_back(bv);
        }
        save_file.data_leg1 = leg1;
        save_file.data_leg2 = leg2;

        write(BWR,"NprFile", save_file);

    }
}


NAMESPACE_END(Grid);
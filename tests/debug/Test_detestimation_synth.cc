/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/debug/Test_detestimation_synth.cc

Copyright (C) 2017

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>
Author: David Murphy <dmurphy@phys.columbia.edu>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

#include <Grid/Grid.h>

using namespace Grid;
using namespace std;

const int max_iter = 5000;
const RealD stop_tol = 1.0e-12;

RealD mean(const std::vector<RealD>& data)
{
    int N = data.size();
    RealD mean(0.0);
    for(int i=0; i<N; ++i){ mean += data[i]; }
    return mean/RealD(N);
}

RealD jack_mean(const std::vector<RealD>& data, int sample)
{
    int N = data.size();
    RealD mean(0.0);
    for(int i=0; i<N; ++i){ if(i != sample){ mean += data[i]; } }
    return mean/RealD(N-1);
}

RealD jack_std(const std::vector<RealD>& jacks, RealD mean)
{
    int N = jacks.size();
    RealD std(0.0);
    for(int i=0; i<N; ++i){ std += std::pow(jacks[i]-mean, 2.0); }
    return std::sqrt(RealD(N-1)/RealD(N)*std);
}

std::vector<RealD> jack_stats(const std::vector<RealD>& data)
{
    int N = data.size();
    std::vector<RealD> jack_samples(N);
    std::vector<RealD> jack_stats(2);

    jack_stats[0] = mean(data);
    for(int i=0; i<N; i++){ jack_samples[i] = jack_mean(data,i); }
    jack_stats[1] = jack_std(jack_samples, jack_stats[0]);
    return jack_stats;
}

template<class Field> class DiagonalOperator  : public LinearOperatorBase<Field> {
public:
  LatticeComplex scale;
  LatticeReal noise;
  DiagonalOperator(GridBase *grid)    : scale(grid), noise(grid)
  {
    GridParallelRNG  pRNG(grid);  
    std::vector<int> seeds({5,6,7,8});
    pRNG.SeedFixedIntegers(seeds);

    random(pRNG,noise);
    RealD amp = 5.0;
    noise =  amp * noise;

    std::vector<Complex> noise_uv(grid->lSites(),1.0);
    unvectorizeToLexOrdArray(noise_uv,noise);

    std::vector<Complex> scale_uv(grid->lSites(),1.0);
  
    
    for (int i = 0; i < grid->lSites()/500; i++) {
      int j = (i+1);
      scale_uv[i]=(4.0*j*j)/(4.0*j*j-1); //= pi/2
      
    }
    for (int i = 0; i < grid->lSites()/500; i+=2) {
      scale_uv[i] = scale_uv[i] * noise_uv[i];
      scale_uv[i+1] = scale_uv[i+1] / noise_uv[i];
    }
    
    vectorizeFromLexOrdArray(scale_uv,scale);


   

  }
  auto diag_product() {
      LatticeComplex logscale(log(scale));
      return exp(sum(logscale));
  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){};
  void OpDiag (const Field &in, Field &out) {};
  void OpDir  (const Field &in, Field &out,int dir,int disp){};

  void Op     (const Field &in, Field &out){
    out = scale * in;
  }
  void AdjOp  (const Field &in, Field &out){
    out = scale * in;
  }
  void HermOp(const Field &in, Field &out){
    double n1, n2;
    HermOpAndNorm(in,out,n1,n2);
  }
  void HermOpAndNorm(const Field &in, Field &out,double &n1,double &n2){
    ComplexD dot;

    out = scale * in;

    dot= innerProduct(in,out);
    n1=real(dot);

    dot = innerProduct(out,out);
    n2=real(dot);
  }
};

template<class Field> class DumbOperator  : public LinearOperatorBase<Field> {
public:
  LatticeComplex scale;
  
  DumbOperator(GridBase *grid, int num_defects)    : scale(grid)
  {
    std::vector<RealD> scale_uv(grid->lSites(),1.0);
    for (int i = 0; i < num_defects; i++) {
    scale_uv[i]=sqrt(2.0);
      
    }
    vectorizeFromLexOrdArray(scale_uv,scale);

  }

  auto diag_product() {
      LatticeComplex logscale(log(scale));
      return exp(sum(logscale));
  }
  auto trace_log() {
      LatticeComplex logscale(log(scale));
      return sum(logscale);
  }
  void OpDirAll  (const Field &in, std::vector<Field> &out){};
  void OpDiag (const Field &in, Field &out) {};
  void OpDir  (const Field &in, Field &out,int dir,int disp){};

  void Op     (const Field &in, Field &out){
    out = scale * in;
  }
  void AdjOp  (const Field &in, Field &out){
    out = scale * in;
  }
  void HermOp(const Field &in, Field &out){
    double n1, n2;
    HermOpAndNorm(in,out,n1,n2);
  }
  void HermOpAndNorm(const Field &in, Field &out,double &n1,double &n2){
    ComplexD dot;

    out = scale * in;

    dot= innerProduct(in,out);
    n1=real(dot);

    dot = innerProduct(out,out);
    n2=real(dot);
  }
};

template<class Field>
void determinant_estimate(GridParallelRNG &RNG, GridCartesian* FGrid, int num_defects, int Nhits) {

    RealD scale = std::sqrt(0.5);
    std::vector<RealD> logdet_D(Nhits);


    ConjugateGradient<Field> CG(stop_tol,max_iter);
    DumbOperator<Field> Mat(FGrid,num_defects);
    
    Field eta(FGrid);
    Field chi(FGrid);
        
    for (int hit = 0; hit < Nhits; hit++) {

        
        gaussian(RNG,eta);
        eta = scale*eta;

        CG(Mat,eta,chi);

        logdet_D[hit] = real(innerProduct(eta,chi)) - norm2(eta);

    }

    std::vector<RealD> logdet_result_vec = jack_stats(logdet_D);

    std::vector<RealD> det_result_vec{logdet_D};
    for (auto&& v : det_result_vec) {
        v = exp(-v);
    }

    std::vector<RealD> det_result = jack_stats(det_result_vec);
    
    size_t power = GridTypeMapper<typename Field::vector_object>::count;
    std::cout << "Number of defects: " << num_defects << " Determinant estimate: " << det_result[0] << " +/- " << det_result[1] << "  Exact result: " << pow(real(Mat.diag_product()()()()),power) << std::endl;

};

int main(int argc, char** argv)
{
    Grid_init(&argc,&argv);
    int threads = GridThread::GetThreads();
    std::cout<<GridLogMessage << "Grid is setup to use "<<threads<<" threads"<<std::endl;

    const int Ls=2;
    const int Nt=4;
    auto latt = GridDefaultLatt();
    latt[3] = Nt;
    //init spacetime grid
    std::cout << GridLogMessage << "Lattice dimensions: "
        << latt << "    Ls: " << Ls << std::endl;

    GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(latt,GridDefaultSimd(Nd,vComplex::Nsimd()),GridDefaultMpi());
    GridCartesian* FGrid = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);


    //init rngs
    
    std::vector<int> seeds5({5,6,7,8});
    GridParallelRNG RNG5(FGrid);
    RNG5.SeedFixedIntegers(seeds5);

    determinant_estimate<LatticeComplex>(RNG5, FGrid,12,30000);
    Grid_finalize();
}

    

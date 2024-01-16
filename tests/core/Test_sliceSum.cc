#include <Grid/Grid.h>



int main (int argc, char ** argv) {
    
    using namespace Grid;

    Grid_init(&argc,&argv);
    GridLogLayout();

    auto latt_size = GridDefaultLatt();
    auto simd_layout = GridDefaultSimd(Nd, vComplex::Nsimd());
    auto mpi_layout = GridDefaultMpi();
    GridCartesian Grid(latt_size, simd_layout, mpi_layout);

    std::vector<int> seeds({1, 2, 3, 4});
    // GridSerialRNG sRNG;
    GridParallelRNG pRNG(&Grid);
    pRNG.SeedFixedIntegers(seeds);

    LatticeComplexD test_data(&Grid);
    gaussian(pRNG,test_data);

    std::vector<TComplex> reduction_reference;
    std::vector<TComplex> reduction_result;

    for (int i = 0; i < Nd; i++) {
        RealD t=-usecond();
        sliceSum(test_data,reduction_reference,i);
        t+=usecond();
        std::cout << " sliceSum took "<<t<<" usecs"<<std::endl;
        
        RealD tgpu=-usecond();
        sliceSumGpu(test_data,reduction_result,i);
        tgpu+=usecond();
        std::cout <<" sliceSumGpu took "<<tgpu<<" usecs"<<std::endl;

    for(int t=0;t<reduction_reference.size();t++){

      std::cout << t<<" reference "<< reduction_reference[t] <<std::endl;
      //TComplex diff = reduced_ref[t]-reduced_gpu[t];
      //assert(abs(TensorRemove(diff)) < 1e-8 );
    }

    for (int t = 0; t<reduction_result.size(); t++) {
        std::cout <<t<<" result " << reduction_result[t] <<std::endl;
    }
    }
    Grid_finalize();
    return 0;
}

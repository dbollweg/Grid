#include <Grid/Grid.h>
#include <Grid/qcd/utils/NprUtils.h>

using namespace std;
using namespace Grid;

int main (int argc, char ** argv) {
  Grid_init(&argc,&argv);
  

  std::string config_prefix;
  if( argc > 1 && argv[1][0] != '-' )
  {
    std::cout<<GridLogMessage <<"Loading from config_prefix "<<argv[1]<<std::endl;

    config_prefix=argv[1];
  }

  std::vector<std::pair<Coordinate, Coordinate> > momenta;
 
  for (int i = 2; i < 3; i++) {
    Coordinate mom1({-i,0,i,0});
    Coordinate mom2({0,i,i,0});
   
    std::pair<Coordinate, Coordinate> mompair(mom1,mom2);
    momenta.push_back(mompair);
  }
  std::vector<int> conf_nums = {100, 1000, 1010};


    for (auto mom: momenta) {   
      std::vector<std::string> conf_names;
      for (auto cnfs: conf_nums) {  
        
    
        std::string filename = config_prefix + "." + std::to_string(cnfs) + "_npr_res_p1_";
        for (auto m1: mom.first){
          filename += std::to_string(m1);
        }

        filename += "_p2_";
        for (auto m2: mom.second){
          filename += std::to_string(m2);
        }

        filename += ".dat";
        conf_names.push_back(filename);
      }
    std::cout << GridLogMessage << "Analyzing " << conf_names.back() << std::endl;
    NPR_analyze npr_reader(conf_names,GridDefaultLatt());
    npr_reader.Test();
    npr_reader.Average();
    }

    return 0;
}
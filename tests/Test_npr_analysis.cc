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
 
  for (int i = 2; i < 4; i++) {
    Coordinate mom1({-i,0,i,0});
    Coordinate mom2({0,i,i,0});
   
    std::pair<Coordinate, Coordinate> mompair(mom1,mom2);
    momenta.push_back(mompair);
  }
  // std::vector<int> conf_nums = {100, 1000, 1010, 1020, 1030, 1040, 1050, 1060, 1070, 1080, 1090};
  // std::vector<int> conf_nums = {765,785,805,825,845};
  std::vector<int> conf_nums = {920, 940, 960, 980, 1000};
  

    for (auto mom: momenta) {   
      std::vector<std::string> conf_names_ml;
      std::vector<std::string> conf_names_ms;
      std::vector<std::string> conf_names_mc;
      for (auto cnfs: conf_nums) {  
        
    
        std::string filename_ml = config_prefix + "." + std::to_string(cnfs) + "_npr_res_ml_p1_";
        std::string filename_ms = config_prefix + "." + std::to_string(cnfs) + "_npr_res_ms_p1_";
        std::string filename_mc = config_prefix + "." + std::to_string(cnfs) + "_npr_res_mc_p1_";
        for (auto m1: mom.first){
          filename_ml += std::to_string(m1);

          filename_ms += std::to_string(m1);
          filename_mc += std::to_string(m1);
        }

        filename_ml += "_p2_";
        filename_ms += "_p2_";
        filename_mc += "_p2_";
        for (auto m2: mom.second){
          filename_ml += std::to_string(m2);
          filename_ms += std::to_string(m2);
          filename_mc += std::to_string(m2);
        }

        filename_ml += ".dat";
        filename_ms += ".dat";
        filename_mc += ".dat";
        conf_names_ml.push_back(filename_ml); 
        conf_names_ms.push_back(filename_ms);
        conf_names_mc.push_back(filename_mc);
      
      }
    std::cout << GridLogMessage << "Analyzing " << conf_names_ml.back() << std::endl;
    NPR_analyze npr_reader_ml(conf_names_ml,GridDefaultLatt());
    // npr_reader.Test();
    npr_reader_ml.Average();
    std::cout << GridLogMessage << "Analyzing " << conf_names_ms.back() << std::endl;
    NPR_analyze npr_reader_ms(conf_names_ms,GridDefaultLatt());
    // npr_reader.Test();
    npr_reader_ms.Average();
    std::cout << GridLogMessage << "Analyzing " << conf_names_mc.back() << std::endl;
    NPR_analyze npr_reader_mc(conf_names_mc,GridDefaultLatt());
    // npr_reader.Test();
    npr_reader_mc.Average();
    }

    return 0;
}
#include <valarray>

class processing_file
{
  int num_feature_elements;
  int num_output_elements;
  int num_patterns;
  int current_pattern_number;
  std::valarray<double> current_feature_vector;
  std::valarray<double> current_target_vector;
  std::valarray<double> mean_feature_vector;
  std::valarray<double> mean_target_vector;
  std::valarray<double> std_feature_vector;
  std::valarray<double> std_target_vector;
  public: 
       read_next_pattern(std::valarray<double> &,std::valarray<double> &);
       processing_file();
       ~processing_file();
};

  

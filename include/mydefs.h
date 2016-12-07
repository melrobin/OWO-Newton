#pragma once
#include <string>
/* ! This struct contains key information about the training, reporting
 * and validation of the MLP.  The reporting piece will likely be strengthened
* greatly in future versions */
typedef struct config_
        {
/* ! Here is the name of the training file */
            std::string trainingfile;
/* ! The number of inputs */
            int inputs;
/* ! The number of outputs */
            int outputs;
/* ! The number of hidden units */
            int hidden_units;
/* ! The number of iterations that the training session will perform */
            int iterations;
/* ! This is currently a misnomer and should be changed to a boolean 
 * variable that states whether or not we want validation */
            int optimal_learning_factor;
/* ! This is the report file string */
            std::string report_file;
/* ! This is the weight file which will be in NuMap format.  We will 
 * document the NuMap format with Doxygen later */
            std::string weight_file;
/* ! For k-fold validation we allow the specification of an integer */
            int k;
            std::string chart_file;
/* ! The name of the validation file */
            std::string val_file;
        } CONFIG ;


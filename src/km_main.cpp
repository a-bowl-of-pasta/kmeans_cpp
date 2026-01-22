
#include "./km_model/kmeans_driver.h"
#include <vector>
#include <iostream>


struct model_config
{

    

    int k_val = default_Kvals[indx_used];     // argv[2]
    int iterations = 100;                    // argv[3]
    double convergence = 0.001;             // argv[4]
    int num_of_runs = 5;                 // argv[5] 



    model_config(int k, int iter, double conv, int totalRuns)
    {   
        k_val = k; 
        iterations = iter; 
        convergence = conv; 
        num_of_runs = totalRuns; 
    }
};

// ================================ test with hard coded values
void inlineRuns()
{

    // hard coded file names for argv[1]
    int indx_used = 3; 
    std::string current_version_out = "v8.0_outputs";

    std::string default_files[] = {"ecoli", "glass", "ionosphere", "iris_bezdek", "landsat", "letter_recognition", "segmentation", "vehicle", "wine", "yeast"}; 
    int default_Kvals[] = {8,6,2,3,6,26,7,4,3,10};
    
    // manually typing file name for argv[1]
    //std::string file_name = "testingMM";           

    // ===== arg variables 

    // so that you can copy and paste :: default_files[indx_used]    file_name
    std::string path_in = "../uci_ml_datasets/with_ground_truth/" + default_files[indx_used] + ".txt";
    std::string path_out = "../outputs/"+ current_version_out +"/"+ default_files[indx_used] + "_out.txt";

}

// =================================== run using args
void actualRun(int argc, char* argv[])
{
     // ===== arg variables 
    std::string file_name; //argv[1]
    int k_val;            // argv[2]
    int iterations;      // argv[3]
    double convergence; // argv[4]
    int num_of_runs;   // argv[5] 


    // ======================= check that args are correct & create k_means API thingy
    try
    { 
        file_name = argv[1]; 
    }
    catch(std::exception e)
    {
        std::cout << "ERROR :: invalid argument | file_Name should be <String>" << std::endl;
        std::cout << "INFO  :: fileName = " << file_name << std::endl;  
        exit(EXIT_FAILURE);        
    }


}


// ============================= main 
int main(int argc, char* argv[])
{

    inlineRuns(); 

return 0; 
}
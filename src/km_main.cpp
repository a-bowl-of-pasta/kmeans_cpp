
#include "./km_model/kmeans_driver.h"
#include <vector>
#include <iostream>


struct model_config
{

    int k_val;     // argv[2]
    int iterations;                    // argv[3]
    double convergence;             // argv[4]
    int num_of_runs;                 // argv[5] 


    model_config(int k, int iter, double conv, int totalRuns)
    {   
        k_val = k; 
        iterations = iter; 
        convergence = conv; 
        num_of_runs = totalRuns;    
    }
    
};



void run_inline() 
{
    // get dataset file location & init model configs 
    int file_array_index = 4; 
    std::string current_version_output = "v9.0_outputs"; 
    
    std::string default_files[10] = {"ecoli", "glass", "ionosphere", "iris_bezdek", "landsat", "letter_recognition", "segmentation", "vehicle", "wine", "yeast"}; 
    int default_Kvals[10] = {8,6,2,3,6,26,7,4,3,10}; 

    std::string path_in = "../uci_ml_datasets/with_ground_truth/" + default_files[file_array_index] + ".txt";
    std::string path_out = "../benchmarks/"+ current_version_output +"/"+ default_files[file_array_index] + "_out.txt";

    model_config configs(default_Kvals[file_array_index], 100, .001, 100); 

    // build model
    km_driver test_model(configs.k_val, configs.num_of_runs, configs.iterations, configs.convergence);
    test_model.build_model(path_in, true, false); 

    // peek dataset data
    std::cout << "\npress <enter> to view \"peek_model_data\" method" << std::endl; 
    std::cin.get(); 
    system("clear"); 

    int num_of_dataset_lines_to_view = 10;
    test_model.peek_model_data(num_of_dataset_lines_to_view);

    // peek data assigned to clusters
    std::cout << "\npress <enter> to view \"peek_cluster_data\" method" << std::endl; 
    std::cin.get(); 
    system("clear"); 


    test_model.peek_cluster_data(configs.k_val);

    // model evaluation summary  
    std::cout << "\npress <enter> to view \"evaluation_summary\" method" << std::endl; 
    std::cin.get();
    system("clear"); 


    test_model.evaluation_summary(); 
}
// ============================= main 
int main(int argc, char* argv[])
{

    run_inline(); 

    return 0; 
}
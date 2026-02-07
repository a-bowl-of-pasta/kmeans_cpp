
#include "./km_model/kmeans_driver.h"
#include <vector>
#include <iostream>
#include <string>


struct model_config
{

    int k_val;              // argv[2]
    int iterations;         // argv[3]
    double convergence;     // argv[4]
    int num_of_runs;        // argv[5] 


    model_config(int k, int iter, double conv, int totalRuns)
    {   
        k_val = k; 
        iterations = iter; 
        convergence = conv; 
        num_of_runs = totalRuns;    
    }
    
};

// ================================= default tester method | test kmeans using some default datsets and configurations 
void test_kmeans(int dataset_choice=3) 
{
    // get dataset file location & init model configs 
    int file_array_index = dataset_choice; 
    std::string current_version_output = "v9.0_outpu ts"; 
    
    std::string default_files[10] = {"ecoli", "glass", "ionosphere", "iris_bezdek", "landsat", "letter_recognition", "segmentation", "vehicle", "wine", "yeast"}; 
    int default_Kvals[10] = {8,6,2,3,6,26,7,4,3,10}; 

    std::string path_in = "../test_datasets/txt_datasets/" + default_files[file_array_index] + ".txt";
    std::string path_out = "../benchmarks/"+ current_version_output +"/"+ default_files[file_array_index] + "_out.txt";
    std::string fileType = "txt"; 

    model_config configs(default_Kvals[file_array_index], 100, .001, 100); 

    // build model
    km_driver test_model(configs.k_val, configs.num_of_runs, configs.iterations, configs.convergence);
    test_model.build_model(path_in, fileType, true, true, false); 

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


    test_model.full_evaluation_summary(); 
}
// ============================= usr sandbox methods

// write some code here
void usr_sandbox(std::string& file_type, int kval, int iterations, double threshold, int totalRuns)
{
    std::string filePath = "../usr_datasets/csv_datasets/penguins.csv";

    model_config config(kval, iterations, threshold, totalRuns);

    km_driver model(config.k_val, config.num_of_runs, config.iterations, config.convergence);
    
    model.build_model(filePath, file_type, true, true);    

    model.internal_evaluation_summary(); 
 
}

// ============================= main 
int main(int argc, char* argv[])
{
    // preview kmeans with already loaded datasets | test_kmeans takes int as param < 0 - 9 >
    test_kmeans(4); 

    // calls user sandbox code
    /*
    std::string file_type; 
    std::cout << "type dataset file type < csv | txt > default is csv :: "; 
    std::cin >> file_type; 

    usr_sandbox(file_type, 5, 100, .001, 100);
    */
    return 0; 
}
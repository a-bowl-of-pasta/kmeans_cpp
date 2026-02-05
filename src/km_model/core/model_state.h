#ifndef MDL
#define MDL

#include <vector>
#include "../cluster.h"
 
struct Model_State
{
    // ===========  model configurations
    int total_runs; 
    double convergence_threshold; 
    int k_value; 
    int max_iterations; 
   
    // =========== per-run tracking | may change per run
    int current_run; 
    int best_run_indx; 
    double best_run_sse; 
    int total_trimmed_dataPoints; 
    int total_padded_dataPoints;
    int total_encrypted_dataPoints; 


    // =========== final evaluation metrics
    double calinski_harabasz_score; 
    double silhouette_score; 
    double jaccard_score; 
    double rand_indx_score;
    
    // =========== data set & ground truth
    std::vector<Data_Point> data_set;  
    std::vector<int> ground_truth_labels; 
 
    // =========== clusters
    std::vector<Cluster> best_run_clust;
    std::vector<Cluster> cluster_list;

    // =========== utility | used in some calculations 
    std::vector<double> max_data_vector; 
    std::vector<double> min_data_vector; 



    // - - - - - - validates and assigns the main values for the constructor 
    void validateData(int k, int iter, double converg, int runs)
    {
        // ----- validates and assigns variables 
        if (k <= 1 || iter <= 0 || runs <= 0) 
        {
           std::cout << "ERROR :: Invalid arguments | values entered : k = " << k 
                     << ", iterations = " << iter
                     << ", total runs ="  << runs <<std::endl; 
           std::cout << "INFO  :: must be : k > 1, iterations > 0, runs > 0 " <<std::endl; 
            exit(EXIT_FAILURE);
        }

        k_value = k; 
        max_iterations = iter; 
        convergence_threshold = converg;
        total_runs = runs; 
    }
   
    Model_State(int k_val, int total_runs, int max_iter, double converg)
    {
        validateData(k_val, max_iter, converg, total_runs);
        
        current_run = 0; 
        best_run_indx = 0; 
        best_run_sse = std::numeric_limits<double>::max(); 
        
        calinski_harabasz_score= std::numeric_limits<double>::min(); 
        silhouette_score = std::numeric_limits<double>::min(); 
        jaccard_score = std::numeric_limits<double>::min();
        rand_indx_score = std::numeric_limits<double>::min();
        total_encrypted_dataPoints = 0; 
        total_padded_dataPoints = 0; 
        total_trimmed_dataPoints = 0; 

        data_set = std::vector<Data_Point>{};   
        ground_truth_labels = std::vector<int>{}; 
        cluster_list = std::vector<Cluster>{};
        best_run_clust = std::vector<Cluster>{};
        max_data_vector = std::vector<double>{};
        min_data_vector = std::vector<double>{};


        cluster_list.reserve(k_val);
        best_run_clust.reserve(k_val);

    }

    Model_State()
    {
        int k = 3; 
        int max_iter = 100; 
        double conv = .001; 
        int runs = 10; 

        validateData(k, max_iter, conv, runs);
        
        current_run = 0; 
        best_run_indx = 0; 
        best_run_sse = std::numeric_limits<double>::max(); 
        
        calinski_harabasz_score = std::numeric_limits<double>::min(); 
        silhouette_score = std::numeric_limits<double>::min(); 
        jaccard_score = std::numeric_limits<double>::min();
        rand_indx_score = std::numeric_limits<double>::min();
        total_encrypted_dataPoints = 0; 
        total_padded_dataPoints = 0; 
        total_trimmed_dataPoints = 0; 

        data_set = std::vector<Data_Point>{};   
        ground_truth_labels = std::vector<int>{}; 
        cluster_list = std::vector<Cluster>{};
        best_run_clust = std::vector<Cluster>{};
        max_data_vector = std::vector<double>{};
        min_data_vector = std::vector<double>{};


        cluster_list.reserve(k);
        best_run_clust.reserve(k);
    }
};

#endif
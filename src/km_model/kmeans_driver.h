#ifndef KM_DRIVER 
#define KM_DRIVER

#include "cluster.h"
#include "./core/model_state.h"
#include "./core/algorithm_backend.h"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <chrono>


class km_io
{

    public:
    

    // - - - - - - basic input / output
    void console_output(std::string output)
    {
        std::cout << output << std::endl; 
    }

    // - - - - - - display the rand index evaluation metric
    void print_rand_index(double rand_index_score)
    {

        if(rand_index_score < 0 || rand_index_score > 1)
        {
            std::cout << "rand index score ::  " << rand_index_score << "\t| ERROR :: rand_index_score ranges 0 to +1, Invalid score :: " << std::endl; 
        }
        else
        {
            if( 0.0 <= rand_index_score && rand_index_score < 0.65)
            {
                std::cout << "rand index score ::  " << rand_index_score << "\t| 0.0 <= rand index score < 0.65, poor ground truth agreement" << std::endl; 
            } 
            else if( 0.65 >= rand_index_score && rand_index_score < 0.80 )
            {
                std::cout << "rand index score ::  " << rand_index_score << "\t| 0.65 <= rand index score < 0.80, weak ground truth agreement" << std::endl; 
            }
            else if(0.80 <= rand_index_score && rand_index_score < 0.90)
            {
                std::cout << "rand index score ::  " << rand_index_score << "\t| 0.80 <= rand index score < 0.90, good ground truth agreement" << std::endl; 
            }
            else if(0.9 <= rand_index_score && rand_index_score <= 1.0)
            {
                std::cout << "rand index score ::  " << rand_index_score << "\t| 0.90 <= rand index score, excellent ground truth agreement" << std::endl; 
            }

        }
        /*
            Ranges from 0 to +1 | closer to +1 is better

            > 0.9: Excellent agreement
            0.8 - 0.9: Good agreement
            0.65 - 0.8: Moderate agreement
            < 0.65: Poor agreement
        */
    }

    // - - - - - - display the jaccard index evaluation metric
    void print_jaccard_index_score(double jaccard_index_score)
    {

        if(jaccard_index_score < 0 || jaccard_index_score > 1)
        {
            std::cout << "jaccard index score ::  " << jaccard_index_score << "\t| ERROR :: jaccard_index_score ranges 0 to +1, Invalid score :: " << std::endl; 
        }
        else
        {
            if( 0.0 <= jaccard_index_score && jaccard_index_score < 0.30)
            {
                std::cout << "jaccard index score ::  " << jaccard_index_score << "\t| 0.0 <= jaccard index score < 0.65, poor ground truth agreement" << std::endl; 
            } 
            else if( 0.30 >= jaccard_index_score && jaccard_index_score < 0.50 )
            {
                std::cout << "jaccard index score ::  " << jaccard_index_score << "\t| 0.65 <= jaccard index score < 0.80, weak ground truth agreement" << std::endl; 
            }
            else if(0.50 <= jaccard_index_score && jaccard_index_score < 0.75)
            {
                std::cout << "jaccard index score ::  " << jaccard_index_score << "\t| 0.80 <= jaccard index score < 0.90, good ground truth agreement" << std::endl; 
            }
            else if(0.75 <= jaccard_index_score && jaccard_index_score <= 1.0)
            {
                std::cout << "jaccard index score ::  " << jaccard_index_score << "\t| 0.90 <= jaccard index score, excellent ground truth agreement" << std::endl; 
            }

        }
        /*
            Ranges from 0 to +1 | closer to +1 is better

            > 0.75: Excellent agreement
            0.5 - 0.75: Good agreement
            0.3 - 0.5: Moderate agreement
            < 0.3: Poor agreement
        */

    }

    // - - - - - - displays the calinski score 
    void print_calinski_score(double ch)
    {
        std::cout << "Calinski Score ::  " << ch << "\t| ranges 0 -> +infinity; the higher the better" << std::endl;
        // Ranges from 0 to +infinity | the higher the better 
    }

    // - - - - - - silhouette score
    void print_silhouette_score(double S_score)
    {
        if(S_score > 1 || S_score < -1)
        {
            std::cout << "Silhouette Score ::  " << S_score << "\t| ERROR :: S_score ranges -1 to +1, Invalid score :: " << std::endl; 
        }
        else
        {
            if( -1.0 <= S_score && S_score < 0)
            {
                std::cout << "Silhouette Score ::  " << S_score << "\t| S-Score < 0.0, bad clustering" << std::endl; 
            } 
            else if( 0 >= S_score && S_score < 0.25 )
            {
                std::cout << "Silhouette Score ::  " << S_score << "\t| 0.0 <= S-Score < 0.25, poor clustering" << std::endl; 
            }
            else if(0.25 <= S_score && S_score < 0.5)
            {
                std::cout << "Silhouette Score ::  " << S_score << "\t| 0.25 <= S-Score < 0.50, weak clustering" << std::endl; 
            }
            else if(0.5 <= S_score && S_score < 0.7)
            {
                std::cout << "Silhouette Score ::  " << S_score << "\t| 0.50 <= S-Score < 0.70, reasonable clustering" << std::endl; 
            }
            else if( 0.7 <= S_score && S_score <= 1)
            {
                std::cout << "Silhouette Score ::  " << S_score << "\t| 0.70 <= S-Score <= 1.0, strong clustering" << std::endl; 
            }
        }
        /* 
            Ranges from -1 to +1 | closer to +1 is better

            > 0.7: Strong clustering
            0.5 - 0.7: Reasonable clustering
            0.25 - 0.5: Weak clustering, some overlap
            < 0.25: No substantial structure / poor clustering
            < 0: Points likely in wrong clusters
        */
    }

    // - - - - - - shows some data assigned to each cluster | help validate data is actually being assigned
   void print_cluster_data(Model_State& state, int clustersToPeek = 0)
   {
        // --- prints min 1 cluster & max all clusters
        if(clustersToPeek <= 0)
        {
            clustersToPeek = 1; 
        }
        else if(clustersToPeek > state.k_value)
        {
            clustersToPeek = state.k_value; 
        }
        std::cout << "total existing clusters :: " << state.cluster_list.size() << std::endl <<std::endl; 

        // --- loop over clusters to peek
        for(int currClust = 0; currClust < clustersToPeek; currClust++ )
        {
            std::vector<int>& clusterData = state.cluster_list.at(currClust).getAssignedData_ref(); 
            
            std::cout << "Cluster < " << state.cluster_list.at(currClust).getID_ref() 
                      << " > first 5 data points | total points in cluster :: "  << clusterData.size() << std::endl; 

            std::cout << std::setfill('-') << std::setw(70)<< ""  << std::setfill(' ') 
                      << std::fixed << std::setprecision(4) << std::endl; 

            // --- print 5 data points per cluster            
            for(int i =0; i < 5; i++)
            {
                Data_Point& point_temp = state.data_set.at(clusterData.at(i)); 
                
                // --- print individual features from feature vector
                for(int feature =0; feature < point_temp.getFeatureDimensions(); feature++)
                {
                    if(feature == 0)
                    {
                        std::cout << point_temp.getFeatureVector_ref().at(feature); 
                    }
                    else
                    {
                        std::cout << std::setw(10) << point_temp.getFeatureVector_ref().at(feature); 
                    }
                }   
                std::cout << std::endl; 
            }

            std::cout  << std::endl; 

        }

   }
   
   // - - - - - - shows a view of the data set | usr defines lines to view, default is whole dataset
    void print_dataset(Model_State& state, int linesToPrint=0)
    {   
        // default, print whole dataset
        if(linesToPrint <= 0)
        { 
            linesToPrint = state.data_set.size();   
        }

        std::cout << "data features"<< std::setw(40) << std::right <<"target clusters" << std::endl; 
        std::cout << std::setfill('-') << std::setw(50)<< ""  << std::setfill(' ') 
                  << std::fixed <<  std::setprecision(4) << std::endl; 

        for(int currPoint = 0; currPoint < linesToPrint; currPoint++)
        {
            Data_Point& point_temp = state.data_set.at(currPoint); 

            for(int feature =0; feature < point_temp.getFeatureDimensions(); feature++)
            {
                if(feature == 0)
                {
                    std::cout << point_temp.getFeatureVector_ref().at(feature); 
                }
                else
                {
                    std::cout << std::setw(8) << point_temp.getFeatureVector_ref().at(feature); 
                }
            }   
            std::cout << std::right << std::setw(20) <<  state.ground_truth_labels.at(currPoint) << std::endl; 
        }

    }

    km_io(){}
};



class km_driver
{

    Alg_Backend backend; 
    Model_State current_state; 
    km_io model_output; 
    
    // - - - - - - finds and sets cluster's mean centroids
    void setNextCentroid()
    {
        // ------ reset clusters with new centroid
        for(int i =0; i < current_state.k_value; i++)
        {   
            
            Cluster& currClust = current_state.cluster_list.at(i); 
            
            // store the mean data vector | call clust::genMeanFeatVector
            std::vector<double> meanCentroid_temp = currClust.mean_feature_vector_classLevel(current_state.data_set); 
            std::string id_temp = "Mean Centroid " + std::to_string(i); 

            Data_Point newPoint_temp(meanCentroid_temp, id_temp); 

            // reset cluster and assign next centroid 
            currClust.resetCluster(); 

            currClust.setCentroid(newPoint_temp); 
        }
    }  
  
    // - - - - - - - starts the next run of the algorithm 
    void start_run(bool peekRuns)
    {
        int current_iteration = 0; 
        bool hasConverged = false; 

        std::vector<double> all_iteration_sse; 

        while (current_iteration < current_state.max_iterations && !hasConverged)
        {
            setNextCentroid(); 
            backend.init_clust(current_state); 

            double iterationSSE = backend.calc_iteration_SSE(current_state);
            all_iteration_sse.push_back(iterationSSE);
            
            // ---- convergence check 
            hasConverged = backend.check_convergance(current_state ,all_iteration_sse, current_iteration, iterationSSE);
            current_iteration++; 
        }        
    }
    
    // - - - - - - - initially sets model & resets between runs
    void set_model(std::string filePath="", bool normalize=false )
    {
            
            if(current_state.current_run == 0)
            {
                for(int i =0; i < current_state.k_value; i++)
                {
                    Cluster clust_temp; 
                    clust_temp.setID(i);  
                    current_state.cluster_list.push_back(std::move(clust_temp));
                }

                backend.load_dataSet(filePath, current_state);
            }
            //~ ------ initialization strategy goes here 
        
            backend.init_forging(current_state); 
        
            // ------ currentRun = 0 is initial run 
            //        decide to normalize or not then increase currentRun from 0 to 1
        
            if(current_state.current_run == 0)
            {
                if(normalize == true)
                    { backend.calc_normalized_data(current_state); }
            
                current_state.current_run++;    
            }
    }
    
    // - - - - - extracts partition from current clustering
    std::vector<int> final_dataset_labels()
    {
        std::vector<int> finalDatasetClustering(current_state.data_set.size(), -1);
        
        // ------ go through each cluster
        for(int label = 0; label < current_state.k_value; label++)
        {
            Cluster& currClust = current_state.cluster_list[label];
            std::vector<int>& clustIndicies = currClust.getAssignedData_ref();
            
            // -------- assign cluster label to each point in current cluster
            for(int i = 0; i < clustIndicies.size(); i++)
            {
                int currFeatureIndx = clustIndicies.at(i);
                finalDatasetClustering[currFeatureIndx] = label;
            }
        }
        
        return finalDatasetClustering;
    }


public: 
    
    // - - - - - - builds the model to be used
    void build_model(std::string& filePath, bool normalize = false, bool printRuns =false)
    {
        set_model(filePath, normalize); 

        //& - - - - - - - - start of timer & main algorithm 
        auto start = std::chrono::high_resolution_clock::now();
        while (current_state.current_run <= current_state.total_runs) // Current run starts at 1 end at total_runs
        {    
            start_run(printRuns);

            set_model();
            
            current_state.current_run++; 
        } 
        auto end = std::chrono::high_resolution_clock::now();
        //& - - - - - - - - end of timer and main algorithm 

        std::chrono::duration<double> seconds_taken = end - start;                
        std::cout << "total time taken for algorithm :: " << seconds_taken.count() << "s" << std::endl;

        current_state.cluster_list = current_state.best_run_clust; 
    }

    // - - - - - - calculate and display the shilouette_score
    void shilouette_score()
    {
        // if there is not a shilouette score, calculate it
        if(current_state.silhouette_score == std::numeric_limits<double>::min())
        {
            current_state.silhouette_score = backend.gen_silhouette_score(current_state); 
        }
        model_output.print_silhouette_score(current_state.silhouette_score); 

    }

    // - - - - - - calculate and display the chi_score
    void calinski_score()
    {   
        // if there is not a calinski score, calculate it
        if(current_state.calinski_harabasz_score == std::numeric_limits<double>::min())
        {
           current_state.calinski_harabasz_score = backend.gen_calinski_index(current_state); 
        }
            model_output.print_calinski_score(current_state.calinski_harabasz_score);  

    }
    
    // - - - - - - calculate and display the rand_index
    void rand_index()
    {
        // if there is not a calinski score, calculate it
        if(current_state.rand_indx_score == std::numeric_limits<double>::min())
        {
            std::vector<int> finalDatasetLabels = final_dataset_labels();
            backend.gen_rand_index_score(current_state, finalDatasetLabels); 
        }

        model_output.print_rand_index(current_state.rand_indx_score); 
        
    }
    
    // - - - - - - calculate and display the jaccard_score
    void jaccard_score()
    {
        // if there is not a calinski score, calculate it
        if(current_state.jaccard_score == std::numeric_limits<double>::min())
        {
            std::vector<int> finalDatasetLabels = final_dataset_labels();
            backend.gen_jaccard_index_score(current_state, finalDatasetLabels); 

        }
        model_output.print_jaccard_index_score(current_state.jaccard_score); 

    }

    // - - - - - - calculate and display the shilouette, chi, rand, and jaccard
    void evaluation_summary()
    {
        // cluster quality metrics | how good are the clusters themselves
        model_output.console_output("\n\nCluster quality evaluation metrics\n"); 
        shilouette_score(); 
        calinski_score(); 

        
        // ground truth comparisons | how closely to gorund truth did the algorithm group the data 
        model_output.console_output("\n\nfinal Kmeans evaluation metrics | ground truth <-> final clustering agreement\n");
        rand_index(); 
        jaccard_score(); 
    }

    // - - - - - - displays data in dataset
    void peek_model_data(int lineToView =0)
    {
        model_output.print_dataset(current_state, lineToView);
    }

    // - - - - - - displays data in assigned to cluster
    void peek_cluster_data(int clustersToView =0)
    {
        model_output.print_cluster_data(current_state, clustersToView); 
    }


    km_driver(int k_val, int total_runs, int max_iter, double converg)
    {
        current_state = Model_State(k_val, total_runs, max_iter, converg);
        backend = Alg_Backend(); 
        model_output = km_io(); 

    }

};



#endif
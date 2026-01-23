
#ifndef ALG_BKND
#define ALG_BKND    

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include "../cluster.h"
#include "./model_state.h"

class Alg_Backend
{

    void create_dataSet_node(const std::string& line, int& id, bool& firstRun, Model_State& state)
    {           
        // read rows element by element
        std::istringstream lineStream(line);
        std::string features; 
        std::vector<double> featureVector_temp;
        std::vector<double>& max_data_vector = state.max_data_vector; 
        std::vector<double>& min_data_vector = state.min_data_vector; 

        int currFeature =0; 
        while(getline(lineStream, features, ' '))
        {
            // ---- skip empty elements 
            if (!features.empty())
            {
                double featToDec = std::stod(features); 
               
                // ----- first run populates min/max vectors, other runs find min/max
                if(firstRun == true)
                {
                    max_data_vector.push_back(featToDec);
                    min_data_vector.push_back(featToDec);
                }
                else
                {
                    // swap attributes if elm in vector is smaller
                    if(max_data_vector.at(currFeature) < featToDec)
                    {
                        max_data_vector.at(currFeature) = featToDec; 
                    }
                    
                    // swap attributes if elm in vector is greater
                    if(min_data_vector.at(currFeature) > featToDec)
                    {
                        min_data_vector.at(currFeature) = featToDec; 
                    }
                }
                featureVector_temp.push_back(featToDec);
                currFeature++; 
            }
        }

        // ---- skip empty vectors  
        if (!featureVector_temp.empty())
        {
            int targetValue = (int) featureVector_temp.at(featureVector_temp.size() - 1);
            featureVector_temp.pop_back(); 

            std::string point_ID_temp = "dataPoint " + std::to_string(id); 
            state.data_set.push_back(Data_Point(featureVector_temp, point_ID_temp));
            state.ground_truth_labels.push_back(targetValue); 
            id++; 

        }              
        
        if(firstRun == true) {firstRun = false; }

    }
 
    bool check_unique(int currClustIndx, int randCentroidIndx, Model_State& state)
    {
        const std::string& dataPointID = state.data_set.at(randCentroidIndx).getPointID(); 
        bool unique = true; 

        for(int i = 0; i < currClustIndx; i++)
        {
            Data_Point& prevCentroids = state.cluster_list.at(i).getCentroid_ref(); 
                 
            if(dataPointID == prevCentroids.getPointID())
            {
                unique = false; 
                break; 
            }
        }

        return unique; 
    }
   
    // - - - - - - gets a mean feature vector using the whole dataset
    std::vector<double> mean_feature_vector_datasetLevel(Model_State& state)
    {
        Data_Point& temp = state.data_set.at(0); 
        int featureDimensions = temp.getFeatureDimensions(); 
    
        std::vector<double> meanFeatureVector(featureDimensions, 0.0); 
    
        int sizeOfData = state.data_set.size();  

        // ------ go through each data vector
        for(int currFeatureVector = 0; currFeatureVector < sizeOfData; currFeatureVector++)
        {
            std::vector<double>& featVector_temp = state.data_set.at(currFeatureVector).getFeatureVector_ref(); 
            
            // -------- calculate the sum of each feature
           for(int feature = 0; feature < featVector_temp.size(); feature++  )
            {
                meanFeatureVector.at(feature) += featVector_temp.at(feature); 
            }

        }

        // ------- divide each feature's sum by total dataPoints
        // ------- this creates a vector of means
        for(int i =0; i < meanFeatureVector.size(); i++)
        {
            meanFeatureVector.at(i) = meanFeatureVector.at(i) / sizeOfData; 
        }

        // ---- return new centroid's dataVector
        return meanFeatureVector; 

    }


    public: 
    
    //  - - - - - - - - Finds the euclidean distance of two points (x1, x2)
    //&                 (X1_0 - X2_0)^2 + (X1_1 - X2_1)^2 + ..... + (X1_n - X2_n)^2
    double sqr_euclid_dist(std::vector<double>& x1_feature_vector, std::vector<double>& x2_feature_vector)
    {
        double finalDistance = 0.0; 
           
        // - - - - loop through x2 features 
        for(int i =0; i < x1_feature_vector.size(); i++)
        {
            // get the features out of the vector 
            double x2_feature = x2_feature_vector.at(i); 
            double x1_feature = x1_feature_vector.at(i); 

            // residual = x1 - x2 
            double curr_sqr_residual = x1_feature - x2_feature; 
           
            // sqr_residual = residual^2
            curr_sqr_residual = curr_sqr_residual * curr_sqr_residual; 

            // finalDist = sum( all sqr_residuals )
            finalDistance += curr_sqr_residual; 
        }

        return finalDistance; 
    }    

    // - - - - - - -  -  assign data from  data_set to cluster 
    void init_clust(Model_State& state) 
    {
        int currentPoint = 0;

        // ----- loops through data set to find closest centroid per dataPoint
        //       Clust saves the dataPoint index from  data_set, not the dataPoint itself
        while(currentPoint < state.data_set.size())
        {
            int bestFitIndex_temp = 0; 
            double currBestFit_temp = std::numeric_limits<double>::max();    

            // - - - finds distance between each centroid and point
            for(int i = 0; i < state.k_value; i ++)
            {
                Data_Point& currCentroid_temp = state.cluster_list.at(i).getCentroid_ref();
                Data_Point& currPoint = state.data_set.at(currentPoint); 
                double distFromCent_temp = sqr_euclid_dist(currPoint.getFeatureVector_ref(), currCentroid_temp.getFeatureVector_ref());
                
                if(distFromCent_temp < currBestFit_temp )
                {
                    bestFitIndex_temp = i; 
                    currBestFit_temp = distFromCent_temp;
                }
               
            }

            state.cluster_list.at(bestFitIndex_temp).assignData(currentPoint); 
            currentPoint++; 
        } 
        
    }    

    // - - - - random data point init
    void init_forging(Model_State& state)
    {
        int currentCluster = 0; 

        // ---- better random numbers
        std::random_device rand_seed; 
        std::mt19937 PRNG(rand_seed());
        std::uniform_int_distribution<int> numRange(0, state.data_set.size() -1); 
        
        // -------- runs until all K clusters have centroids
        while (currentCluster < state.k_value)
        {
            int centroidIndx = numRange(PRNG); 

            // ---- first cluster, no need to compare
            if(currentCluster == 0)
            {
                Cluster& currClust = state.cluster_list.at(currentCluster);                 
                currClust.setCentroid(state.data_set.at(centroidIndx)); 
                currentCluster++; 
            }
            else // ---- compare with other clusters, no repeated centroids 
            {               
                // ---- assign centroid if unique
                if(check_unique( currentCluster, centroidIndx, state))
                {
                    Cluster& currClust = state.cluster_list.at(currentCluster); 
                    currClust.setCentroid(state.data_set.at(centroidIndx));  
                    currentCluster++;
                }
            }
        }

        init_clust(state); 
    }

    // - - - - - - reads input file, 
    void load_dataSet( const std::string&  dataFilePath, Model_State& state)
    {
        std::string line; 
        std::ifstream fn(dataFilePath); 
        bool isFirstLine = true; 
        int id = 0; 
        bool firstRun = true; 

        // ---- check if file is open
        if(!fn.is_open())
        {
            std::cout << "ERROR :: file is not opening"<<std::endl; 
            std::cout << "INFO  :: file path = " << dataFilePath <<std::endl; 
            exit(EXIT_FAILURE);
        }

        // -------- read lines from file
        while(getline(fn, line))
        {
            // ---- first line = header, following lines = data 
            if(isFirstLine == false)
            {
               create_dataSet_node( line, id, firstRun, state); 
            }
            else 
            {
                id++; 
                isFirstLine = false;    
            }
        }
       
        fn.close(); 

    }

    // - - - -  data normalization
    void calc_normalized_data(Model_State& state)
    {

        // ------ goes through each feature vector
        for(int currFeatVectIndx =0; currFeatVectIndx < state.data_set.size(); currFeatVectIndx++)
        {
            std::vector<double>& currFeatVect_temp = state.data_set.at(currFeatVectIndx).getFeatureVector_ref(); 
            
            // ----- goes through each elm in feature vector
            for(int featureIndx = 0; featureIndx < currFeatVect_temp.size(); featureIndx++ )
            {
                double currFeature_temp = currFeatVect_temp.at(featureIndx); 
                double maxVal_temp = state.max_data_vector.at(featureIndx); 
                double minVal_temp = state.min_data_vector.at(featureIndx); 
                
                // x' = x - min(x) / max(x) - min(x)
                double normalizedFeature = 0.0; 
                if((maxVal_temp - minVal_temp) != 0)
                {
                    normalizedFeature = (currFeature_temp - minVal_temp) / (maxVal_temp - minVal_temp);
                }

                currFeatVect_temp.at(featureIndx) = normalizedFeature; 
            }
        }
    }

    // - - - - finds the iteration SSE
    double calc_iteration_SSE(Model_State& state)
    {
        double totalIterationSSE = 0.0;
        for(int i =0; i < state.k_value; i++)
        {
            state.cluster_list.at(i).sum_squared_error(state.data_set); 
            totalIterationSSE += state.cluster_list.at(i).getSSE_classLevel(); 
        }
        return totalIterationSSE; 
    }
   
    // - - - - creates a new node (dataPoint) for the dataSet vector | not used outside this class
    bool check_convergance(Model_State& state, std::vector<double>& iter_sse_vector, int current_iteration, double iterSSE)
    {
        bool converged = false; 
        
        // ---- convergence check for :: i > 0 or i == 0
        if(current_iteration > 0)
        { 
            // (SSE^t-1 - SSE^t) / SSE^t-1
            double prevSSE = iter_sse_vector.at(current_iteration - 1);
            double currentSSE = iter_sse_vector.at(current_iteration); 
            double convCheck = (prevSSE - currentSSE) / prevSSE;

            if (convCheck < state.convergence_threshold)
            {
                converged = true; 
            }
        }
        else 
        {
            if(iterSSE < state.convergence_threshold) 
            {              
                converged =  true; 
            }
        }

        // ---- if converged, compare current iteration to algorithm's best iteration
        if(state.best_run_sse > iterSSE && converged == true)
        {
            // sets values for the current best 'run' 
            state.best_run_sse = iterSSE; 
            state.best_run_indx = state.current_run; 

            // deep copy best run clusters
            state.best_run_clust = state.cluster_list;             
        } 

        return converged; 
    }
   
    // - - - - - - CH = (SSB / (k - 1)) / (SSW / (n - k))
    double gen_calinski_index(Model_State& state)
    {
        // between cluster sum of squares = SSB 
        double between_cluster_sum_of_squares = 0.0; 
        std::vector<double> meanDataVec =  mean_feature_vector_datasetLevel(state); 
       
        for(int i = 0; i < state.k_value; i++)
        {
            Cluster& currClust = state.cluster_list.at(i); 
            Data_Point& currCent = currClust.getCentroid_ref(); 

            // gets the squared euclidean distance | (centroid - meanVec)^2
            double pt1 = sqr_euclid_dist(currCent.getFeatureVector_ref(), meanDataVec); 

            between_cluster_sum_of_squares += currClust.cluster_size() * pt1; 
        }

        double pt1 = between_cluster_sum_of_squares / (state.k_value - 1); 
        double pt2 = state.best_run_sse / (state.data_set.size() - state.k_value); 

        return pt1 / pt2; 
    }

    // - - - - - - silhouette score 
    double gen_silhouette_score(Model_State& state)
    {
        double currentClustSilhouette = 0.0 ;

        for(int i = 0; i < state.k_value; i++)
        {    
            Cluster& currClust = state.cluster_list.at(i); 
            
            currentClustSilhouette += currClust.silhouette_score(state.data_set, state.cluster_list);
            
        }
        return (currentClustSilhouette / state.k_value); 

    }

    // - - - - - calculates Rand Index between two partitions | good for balanced datasets - no disproportionate FN,FP,TP,TN
    //&          (truePositives + trueNegatives) / total pairs
    void gen_rand_index_score(Model_State& state, std::vector<int>& datasetLabels)
    {
        int totalPoints = datasetLabels.size();
        int truePositives = 0;  // same final clust in preditction & ground truth
        int trueNegatives = 0;  // different final clust in preditction & ground truth
        
        // ------ compare all pairs
        for(int x1 = 0; x1 < totalPoints; x1++)
        {
            for(int x2 = x1 + 1; x2 < totalPoints; x2++)
            {
                bool sameInPrediction = (datasetLabels[x1] == datasetLabels[x2]);
                bool sameInGroundTruth = (state.ground_truth_labels[x1] == state.ground_truth_labels[x2]);
                
                // correctly predicted two points in the same cluster
                if(sameInGroundTruth && sameInPrediction)
                    { truePositives++;}

                // correctly predicted two points in different clusters
                else if(!sameInGroundTruth && !sameInPrediction)
                    { trueNegatives++; }
            }
        }
        
        double totalPairs = (double)(totalPoints * (totalPoints - 1)) / 2;
        double finalScore = (double)((truePositives + trueNegatives) / totalPairs);

        state.rand_indx_score = finalScore;
    }

    // - - - - - calculates Jaccard Index between two partitions | good for disproportionate datasets - far more TN than the rest
    //&          true positives / (true positives + false negative + false positives)
    void gen_jaccard_index_score(Model_State& state, std::vector<int>& datasetLabels)
    {
        int totalPoints = datasetLabels.size();
        int truePositives = 0;  // same final clust in preditction & ground truth
        int falsePositives = 0;  // same final clust in preditction but different in ground truth
        int flaseNegatives = 0;  // different final clust in preditction but same in  ground truth
        
        // ------ compare all pairs
        for(int x1 = 0; x1 < totalPoints; x1++)
        {
            for(int x2 = x1 + 1; x2 < totalPoints; x2++)
            {
                bool sameInPrediction = (datasetLabels[x1] == datasetLabels[x2]);
                bool sameInGroundTruth = (state.ground_truth_labels[x1] == state.ground_truth_labels[x2]);
                
                // correctly predicted two points in same cluster
                if(sameInPrediction && sameInGroundTruth)
                    { truePositives++; }

                // incorrectly predected two points in the same cluster
                else if(sameInPrediction && !sameInGroundTruth)
                    { falsePositives++; }

                // 
                else if(!sameInPrediction && sameInGroundTruth)
                    { flaseNegatives++;}
            }
        }

        // does not account for true negatives | more strict but better for disproportionate sets
        //~ disproportinate where there may be 80 true positives but only 5 true negatives
        double totalPairs = (double) truePositives + falsePositives + flaseNegatives; 
        double finalScore = (double)(truePositives / totalPairs); 

        state.jaccard_score = finalScore;
    }


    Alg_Backend()
    {}
};

#endif
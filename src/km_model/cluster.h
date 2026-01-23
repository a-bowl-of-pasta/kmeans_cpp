#ifndef CLUST_H
#define CLUST_H

#include <iostream>
#include <string>
#include <vector>
#include <limits>

class Data_Point
{

    // data ID | distance from nearest cluster | data
    std::string dataPoint_class_id; 
    std::vector<double> feature_vector;
   
    public : 
    // = = = = = = = = = = = = = = setters and getters 

    // - - - - - - - - - - data point ID 
    void setPointID(std::string& id) 
        { dataPoint_class_id = std::move(id); }
    
    std::string& getPointID()
        { return dataPoint_class_id;}

    
    // - - - - - - - - - - - data vector 
    void setFeatureVector(std::vector<double>& updatedVector)
        { feature_vector = std::move(updatedVector); }
    
    std::vector<double> getFeatureVector_copy()
        { return feature_vector;}

    
    std::vector<double>& getFeatureVector_ref()
        { return feature_vector; }    

    int getFeatureDimensions()
        { return feature_vector.size(); }


    // - - - -  squared euclidean distance calculation
    //&         (X1_0 - X2_0)^2 + (X1_1 - X2_1)^2 + ..... + (X1_n - X2_n)^2
    double sqr_euclid_dist(Data_Point& x2_point)
    {
        double finalDistance = 0.0; 
        std::vector<double>& x2_feature_vector = x2_point.getFeatureVector_ref(); 
           
        // - - - - loop through x2 features 
        for(int i =0; i < feature_vector.size(); i++)
        {
            // get the features out of the vector 
            double x2_feature = x2_feature_vector.at(i); 
            double x1_feature = feature_vector.at(i); 

            // residual = x1 - x2 
            double curr_sqr_residual = x1_feature - x2_feature; 
           
            // sqr_residual = residual^2
            curr_sqr_residual = curr_sqr_residual * curr_sqr_residual; 

            // finalDist = sum( all sqr_residuals )
            finalDistance += curr_sqr_residual; 
        }

        return finalDistance; 
    }    

    Data_Point(std::vector<double>& feature_vector, std::string& instance_id)
    :   dataPoint_class_id(std::move(instance_id)),
        feature_vector(std::move(feature_vector))
    {}
    
    Data_Point()
    :   dataPoint_class_id(),
        feature_vector()
    {}

};

class Cluster
{

    double class_level_SSE;
    double class_level_SScore; 
    std::string clust_id; 
    Data_Point centroid; 
    std::vector<int> clust_data_indicies;

    public: 

    // - - - - - - - - - - - - indicies of datapoints assigned to cluster
    void assignData(int dataSetIndex)
        { clust_data_indicies.push_back(dataSetIndex);}


    std::vector<int> getAssignedData_copy()
        {return clust_data_indicies; }

    std::vector<int>& getAssignedData_ref()
        {return clust_data_indicies; }


    int cluster_size()
        {return clust_data_indicies.size();}


    // - - - - - - - - - - -  centroid 
    void setCentroid(Data_Point& cent) 
        { centroid = cent;}
    
    Data_Point& getCentroid_ref()
        {return centroid; }

    Data_Point getCentroid_copy()
        {return centroid;}


    // - - - - - - - - - - - cluster ID
    void setID(int id_num)
        { clust_id = "c" + std::to_string(id_num); }

    std::string& getID_ref()
        {return clust_id;}
    
    std::string getID_copy()
        {return clust_id;}

    
    // - - - - - - - - - - - class SSE
    double getSSE_classLevel()
        {return class_level_SSE;}


    // - - - - - - - - - - - class shiloette score 
    double getShilScore_classLevel()
        {return class_level_SScore;}
  
    
    // - - - - - - - - - - - reset clust to default values
    void resetCluster()
    {
        centroid = Data_Point{};
        class_level_SSE = 0.0; 
        class_level_SScore = 0.0; 
        clust_data_indicies.clear();
    }

    // - - - - - - - - -  square euclidean distance 
    //&                   (X1_0 - X2_0)^2 + (X1_1 - X2_1)^2 + ..... + (X1_n - X2_n)^2
    double dist_from_centroid(std::vector<double>& x2_feature_vector)
    {
        double finalDistance = 0.0; 
           
        // - - - - loop through x2 features 
        for(int i =0; i < centroid.getFeatureVector_ref().size(); i++)
        {
            // get the features out of the vector 
            double x2_feature = x2_feature_vector.at(i); 
            double x1_feature = centroid.getFeatureVector_ref().at(i); 

            // residual = x1 - x2 
            double curr_sqr_residual = x1_feature - x2_feature; 
           
            // sqr_residual = residual^2
            curr_sqr_residual = curr_sqr_residual * curr_sqr_residual; 

            // finalDist = sum( all sqr_residuals )
            finalDistance += curr_sqr_residual; 
        }

        return finalDistance; 
    }    


    // - - - - - - - - - - - - - -  - SSE calculation 
    void sum_squared_error( std::vector<Data_Point>& dataSet)
    {
        // - - initialize variables
        int totalPointsInCluster = clust_data_indicies.size(); 

        double clustSSE = 0; 

        // - - goes through each data point assigned to cluster 
        for(int currPoint = 0; currPoint < totalPointsInCluster; currPoint++ )
        {
            int pos = clust_data_indicies.at(currPoint); 
            std::vector<double>& currFeaturesVector = dataSet.at(pos).getFeatureVector_ref(); 

            double currFeatureVectorSSE = dist_from_centroid(currFeaturesVector);

            clustSSE += currFeatureVectorSSE; 
        }

        class_level_SSE = clustSSE; 
    }


    // - - - - - - finds the mean data vector for the data assined to this cluster
    std::vector<double> mean_feature_vector_classLevel( std::vector<Data_Point>& dataSet)
    {
        
        // - - initialize variables | finalMeanVector - vector of each column's mean
        int featureVectorDimensions = dataSet.at(0).getFeatureDimensions(); 
    
        std::vector<double> finalMeanVector(featureVectorDimensions, 0.0); 
    
        // ~ cluster_Data_indicies - vector of indicies for each point assigned to the cluster
        int totalPointsInCluster = clust_data_indicies.size(); 

        // - - goes through data assigned to cluster and sums each column
        for(int currentPoint = 0; currentPoint < totalPointsInCluster; currentPoint++)
        {
            int pos = clust_data_indicies.at(currentPoint); 
    
            std::vector<double>& features_currentPoint = dataSet.at(pos).getFeatureVector_ref(); 
            
            // sums columns 
           for(int currFeature = 0; currFeature < features_currentPoint.size(); currFeature++)
            {
                finalMeanVector.at(currFeature) += features_currentPoint.at(currFeature); 
            }
        }

        // divide each column's sum by the total amount of points
        for(int i =0; i < finalMeanVector.size(); i++)
        {
            finalMeanVector.at(i) = finalMeanVector.at(i) / totalPointsInCluster; 
        }

        // return vector of column means
        return finalMeanVector; 
    }


    // - - - - - - cohesion using a comprehensive point to point strategy
    //&             a(i) = (currPoint - every_Other_Point_In_Clust)^2
    std::vector<double> clust_cohesion_pairwise(std::vector<Data_Point>& dataset)
    {
        int currClustTotalPoints = cluster_size(); 

        std::vector<double> allCohesionScores; 

        // --- loop through curr clust assigned points
        for(int currPoint = 0; currPoint < currClustTotalPoints; currPoint++)
        {
            double currPointCohesionCoeff = 0.0; 
            int x1_pos = clust_data_indicies.at(currPoint); 

            // --- (x_1 - x_2)^2 | x_2 = each point assigned to clust excluding x_1
            for(int allPoints = 0; allPoints < currClustTotalPoints; allPoints++)
            {
                int x2_pos = clust_data_indicies.at(allPoints); 
                if(x1_pos != x2_pos)
                {
                    currPointCohesionCoeff += dataset.at(x1_pos).sqr_euclid_dist(dataset.at(x2_pos));         
                }
            }
            allCohesionScores.push_back((currPointCohesionCoeff / (currClustTotalPoints -1) )); 

        }
        
        return allCohesionScores;
    }


    // - - - - - - separation using point to point, pairwise, strategy
    //&            b(i) = (currPoint - every_Other_Point_In_other_clusters)^2
    void separation_pairwise(std::vector<double>& b_i, std::vector<Data_Point>& dataset, std::vector<int>& clustB_data)
    {
        int clustA_TotalPoints = cluster_size(); 
        int clustB_TotalPoints = clustB_data.size(); 

        // --- loop through clustA assigned points
        for(int currPoint = 0; currPoint < clustA_TotalPoints; currPoint++)
        {
            double pointSeparation = 0.0; 
            int x1_pos = clust_data_indicies.at(currPoint); 

            // --- (x_1 - x_2)^2 | loop through clustB
            for(int allPoints = 0; allPoints < clustB_TotalPoints; allPoints++)
            {
                int x2_pos = clustB_data.at(allPoints); 
        
                pointSeparation += dataset.at(x1_pos).sqr_euclid_dist(dataset.at(x2_pos));         
            }
            pointSeparation = pointSeparation / clustB_TotalPoints; 

            b_i.at(currPoint) = std::min(b_i[currPoint], pointSeparation);
        }        

    }

    
    // - - - - - - - final shilouette score 
    //&              s(i) = b(i) - a(i) / max(a(i), b(i))
    double silhouette_score(std::vector<Data_Point>& data_set, std::vector<Cluster>& clust_list)
    {
        // a(i) = indicidual cohesion && b(i) = individual separation 
        std::vector<double> individual_cohesion = clust_cohesion_pairwise(data_set); 
        std::vector<double> individual_separation(cluster_size(), std::numeric_limits<double>::max()); 

        // ---- loop through k cluster and gets the separation values
        for(int currClust = 0; currClust < clust_list.size(); currClust++)
        {
            // -- skip this cluster
            if(clust_list.at(currClust).getID_ref() != clust_id)
            {
                separation_pairwise( individual_separation,
                                     data_set, 
                                     clust_list.at(currClust).getAssignedData_ref()); 
            }
        }

        //  s(i) = b(i) - a(i) / max(a(i), b(i))
        double finalSilhouette = 0.0; 
        for(int i =0; i < cluster_size(); i++)
        {
            double a_i = individual_cohesion.at(i); 
            double b_i = individual_separation.at(i); 
            double s_i = (b_i - a_i) / std::max(a_i, b_i);

            finalSilhouette += s_i; 
        } 

        class_level_SScore = (finalSilhouette / cluster_size()); 
        return class_level_SScore ; 
    }

    
    Cluster()
    {
        class_level_SSE  = 0.0;
        class_level_SScore = 0.0; 
        clust_id = ""; 
        centroid = Data_Point{}; 
        clust_data_indicies = std::vector<int>{};
    }


};




#endif
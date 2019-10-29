## 1D random spines for-loop:
for(j in 1:5000){                                                        # generate list of randomly generated attachment points to calculate random clustering
  random_spines_1D <- sample(possible_dist, total_spines)              # take sample from possible locations on dendrite to match total number of spines
  test_dist_1D <- as.matrix(dist(random_spines_1D))                    # create distance matrix for random sample
  UPGMA_test_1D <- IdClusters(test_dist_1D,                            # run UPGMA with randomized data
                              method = "UPGMA",                        # used UPGMA method
                              cutoff = 0.75,                           # cutoff from Yadav paper
                              showPlot = FALSE)                        # suppresse dendrogram
  UPGMA_test_1D$rn <- rownames(UPGMA_test_1D)                          # add row name (spine ID) to cluster ID number, not sure if necessary  
  cluster_all_test_1D <- cbind(UPGMA_test_1D, random_spines_1D)        # add cluster number/spine ID to randomly generated soma distances
  cluster_all_test_1D <- cluster_all_test_1D[order(cluster_all_test_1D$cluster), ] # order by cluster number
  cluster_freq_test_1D <- table(cluster_all_test_1D$cluster)           # generate table with number of spines in each cluster  
  cluster_freq_test_1D <- as.data.frame(cluster_freq_test_1D)          # change table to data frame 
  cluster_freq_test_1D$is_clustered <- as.numeric(cluster_freq_test_1D$Freq > 1) # return 1 if  > 1 spine in cluster, 0 if not-- 1s reflect true "cluster"
  num_clusters_test_1D <- sum(cluster_freq_test_1D$is_clustered,                 # count # of clusters on this segment
                              na.rm = TRUE) 
  num_clusters_test_1D[is.na(num_clusters_test_1D)] <- 0               # changes possible NA to 0
  num_clusters_test_1D <- as.numeric(num_clusters_test_1D)             # convert to numerical form
  test_num_clusters_all_1D <- rbind(test_num_clusters_all_1D,          # add total number of clusters to running list
                                    num_clusters_test_1D) 
  spines_clustered_test_1D <- sum(cluster_freq_test_1D$Freq > 1,       # sum total number of spines that are in clusters of > 1 spine
                                  na.rm = TRUE) 
  spines_not_test_1D <- data.frame(total_spines - spines_clustered_test_1D) # calculate how many spines are not clustered
  spines_clustered_test_1D[is.na(spines_clustered_test_1D)] <- 0            # return 0 instead of Na if no spines are clustered in the random sample
  spines_clustered_test_1D <- data.frame(spines_clustered_test_1D)          # convert to data.frame
  test_spines_clustered_all_1D <- rbind(test_spines_clustered_all_1D, 
                                        spines_clustered_test_1D)           # add number of clustered spines to running list 
} # end of 1D random spines for-loop

#Determining 1D Cscore--the closer to 1, the higher probability that a given number has more clustered spines than random normal distribution
test_spines_clustered_all_1D <- as.matrix(test_spines_clustered_all_1D) #changes data.frame to matrix (not sure if this is neccesary but it works)
test_spines_clustered_all_1D <- as.numeric(test_spines_clustered_all_1D) #changes all vales to numeric (again, not sure if neccessary)
std_test_1D <- sd(test_spines_clustered_all_1D) #generates STD of total number of spines clustered in entire random sample
mean_test_1D <- mean(test_spines_clustered_all_1D) #generates average number of spines in a cluster in random data
curve_dnorm_1D <- dnorm(test_spines_clustered_all_1D, mean_test_1D, std_test_1D) #gives probability density function, or height of probability distribution at each point(height = frequency)
std_curve_1D <- sd(curve_dnorm_1D)
mean_curve_1D <- mean(curve_dnorm_1D)
Cscore_1D <- pnorm(spines_clustered_1D, mean_test_1D, std_test_1D) #calculates the probability that given the random sample distribution (curve mean+SD) the total number of spines observed is higher

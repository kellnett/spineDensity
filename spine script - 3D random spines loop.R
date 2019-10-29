# 3D random spines for loop
#  for(j in 1:50){
#    test_data_X <- data.frame(sample(df$X), df$Y, df$Z) # randomize the X's, Y's, and Z's to make a "biologically plausible" dataframe.
#    colnames(test_data_X) <- c( "x", "Y", "Z")
#    test_data_Y <- data.frame(df$X, sample(df$Y), df$Z)
#    colnames(test_data_Y) <- c( "x", "Y", "Z")
#    test_data_Z <- data.frame(df$X, df$Y, sample(df$Z))
#    colnames(test_data_Z) <- c( "x", "Y", "Z")
#    test_data   <- rbind(test_data_X, test_data_Y, test_data_Z)
#    test_data_final <- data.frame(sample_n(test_data, total_spines))
#    test_dist_3D <- as.matrix(dist(test_data_final)) #creates distance matrix for random sample

#    UPGMA_test_3D <- IdClusters(test_dist_3D, method = "UPGMA", cutoff = 0.75, showPlot = FALSE) #run cluster analysis with cutoff used in Yadav paper
#    UPGMA_test_3D$rn <- rownames(UPGMA_test_3D)   #adds a column with rownows to keep spine ID, not sure if this step is neccessary
#    cluster_all_test_3D <- cbind(UPGMA_test_3D, test_data_final)   #get the coordinates for each spine ID, if kept in row name/number order, will correctly correspond to each spine
#   cluster_all_test_3D <- cluster_all_test_3D[order(cluster_all_test_3D$cluster), ] #sort data frame by cluster number
#    cluster_freq_test_3D <- table(cluster_all_test_3D$cluster) #create a table counting how many times each cluster Variable occurs (i.e. how many spines in each cluster)
#    cluster_freq_test_3D <- as.data.frame(cluster_freq_test_3D)
#    cluster_freq_test_3D$is_clustered <- as.numeric(cluster_freq_test_3D$Freq > 1) # create column where 1 means there is more than one spine in a cluster or 0 if just 1
#    num_clusters_test_3D <- sum(cluster_freq_test_3D$is_clustered, na.rm = TRUE)
#    num_clusters_test_3D[is.na(num_clusters_test_3D)] <- 0
#    num_clusters_test_3D <- data.frame(num_clusters_test_3D)
#    test_num_clusters_all_3D <- rbind(test_num_clusters_all_3D, num_clusters_test_3D)
#    spines_clustered_test_3D <- cluster_freq_test_3D %>% 
#      group_by(is_clustered) %>% 
#      filter(is_clustered == 1) %>% 
#      summarise(spines_clustered_test_3D = sum(Freq))
#    ifelse(is.null(spines_clustered_test_3D) == TRUE, spines_clustered_test_3D <- 0, spines_clustered_test_3D <- spines_clustered_test_3D[1,2])
#    spines_not_test_3D <- as.numeric(total_spines - spines_clustered_test_3D)
#    spines_clustered_test_3D[is.na(spines_clustered_test_3D)] <- 0
#    spines_clustered_test_3D <- data.frame(spines_clustered_test_3D)
#    test_spines_clustered_all_3D <- rbind(test_spines_clustered_all_3D, spines_clustered_test_3D)

#  } # end of random spines 3D for-loop


#  test_spines_clustered_all_3D <- as.numeric(as.matrix(test_spines_clustered_all_3D))
#  std_test_3D <- sd(test_spines_clustered_all_3D)
#  mean_test_3D <- mean(test_spines_clustered_all_3D)
#  curve_dnorm_3D <- dnorm(test_spines_clustered_all_3D, mean_test_3D, std_test_3D)
#  std_curve_3D <- sd(curve_dnorm_3D)
#  mean_curve_3D <- mean(curve_dnorm_3D)
#  Cscore_3D <- pnorm(spines_clustered_3D, mean_test_3D, std_test_3D)


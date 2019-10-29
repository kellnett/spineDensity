for(k in 1:5000) {  ## random ac for-loop, use random distances from soma to generate list of randomized ac values
  test_data <- sample(possible_dist, total_spines) # take sample equal to total number of spines from all possible soma distnaces 
  test_dist <- as.matrix(dist(test_data))          # create distance matrix of random soma distances
  assign("test_cluster",                           # assign the following to the output file 'test_cluster':
         agnes(test_data,                        # values from agnes function (package 'cluster') on the random soma distances
               metric = "euclidean",                   # use euclidean distance
               method = "average"))                    # use the UPGMA method
  ac_test   <- rbind(ac_test, test_cluster$ac)     # add row to running list of test ac's
} ## end of random ac for-loop

cScore_ac_1D    <- sum(as.numeric(ac_test < dendrite_ac))/5000 # average (divide by 1000 samples) how many times random ac is smaller than dendrite_ac, the smaller the value, the more "clustering"
df$cScore_ac_1D <- cScore_ac_1D                                # add a row to df with the cScore of the dendrite
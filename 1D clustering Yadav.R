## Code Overview
## I. Intro
## II. Initializing code
## III. FileList for-loop
## IV. Compile and Summarise Data
## V. Save and write .csv files

#  This code models Yadav et al. 2012 analysis of dendritic spines using 
#  unsupervised agglomerative hierarchical clustering method. Initial 
#  creation of this code was by S. Johnson and was tailored towards 
#  clustering of attachment points

#  Recent updates made by K. Nett in collaboration with J. Kent to:
#  Improve calculations/readability and establish clustering method that 
#  uses 3D coordinates (spine heads)

#   Various clustering methods tested: 

#   Method coded by S.Johnson looking at attachment points 
#   (i.e. SOMA.DISTANCE) using 'agnes' function (package 'cluster')
#   Referred to as "dendrite_ac" within the code

#   Method still utilizing attachment points (i.e. SOMA.DISTANCE) 
#   but using function 'IdClusters' from the package 'DECIPHER'
#   Referred to as "Cscore_ac_1D"

#   Method utilizing attachment point 3D coordinates but using
#   function 'IdClusters' from the package 'DECIPHER'
#   Referred to as "Cscore_1D"

#   Method looking at clustering of 3D spine head locations using 
#   function 'IdClusters' from the package 'DECIPHER'
#   Referred to as "Cscore_3D"

#   Nearest neighbor analysis using both 1D attachment points and 
#   3D spine head locations (nn_dist_1D or nn_dist_3D)


#   Package Installation instructions (you only need to install these once)
#   enter install.packages("package name") into console of Rstudio,
#   example: install.packages("dplyr")
#   Special cases: bio3d and DECIPHER
#   bio3D: install.packages("bio3d", dependencies = TRUE)
#   DECIPHER: source("https://bioconductor.org/biocLite.R") 
#             biocLite("DECIPHER")
#             (when prompted, enter 'a' for "all")


timer <- proc.time() # start timing code


# Locate packages in library

library(dplyr)           
library(ggplot2)        
library(lattice)
library(latticeExtra)
library(cluster)
library(stats)
library(readr)
library(tidyr)
library(bio3d)
library(DECIPHER)
library(svDialogs)

# Initialle Code: Select file(s) and define group
      
fileChosen <- file.choose()                           # opens a folder dialog box to select first file in folder containing data from 1 treatment group 
filePath   <- dirname(fileChosen)                     # gets the directory name that that file is in
setwd(filePath)                                       # sets the working directory. important for saving files later
fileList   <- list.files(filePath)                    # list the files in the directory to go through
fileList   <- grep(".csv", fileList, value = TRUE)    # finds only lists .csv files in case there are other types
data_all   <- data.frame()                            # initializes a data.frame that will store all of the data for the experiment
group_name <- dlgInput("Enter Group Name", 
                       Sys.info()["user"])$res        # prompt given in console to enter the animal treatment group (i.e. Sal-D0)

# File list for-loop: Apply code to each file in fileList

for(i in 1:length(fileList))  {                             # peform code inside {} for each .csv file in the file list
  fileName     <- unlist(strsplit(fileList[i], "[.]"))[1]   # remove '.csv' from the file name
  fileLoc      <- file.path(filePath, fileList[i])          # create path to file (will differ from initial FileChosen as loop progresses)
  df           <- read_csv(fileLoc, col_types = cols())     # read the csv file, uses column names from NeuronStudio file
  colnames(df) <- gsub("-", "_", colnames(df))              # replaces dashes ('-') w/ underscores ('_') for column names (R doesn't like dashes)
  df           <- df[complete.cases(df$SOMA_DISTANCE) ,]    # removes cases where there is no soma distance data
  df           <- df[complete.cases(df$RAYBURST_VOLUME), ]  # " " spine volume data
  df           <- df[complete.cases(df$MAX_DTS), ]          # " " length data
  df$file      <- fileName                                       # add a column identifying the file name for each spine (row)
    
## Current NeuronStudio data should only have 1 section/dendrite; however, it is possible
## to have multiple if they are not linked together. The next 6 lines calculate total
## length when there is more than 1 section
    
total_length <- df %>%                                    # dplyr package: define total length of dendritic segment
  group_by(SECTION_NUMBER)  %>%                           # group data frame by section number
  summarise(section_length = max(SECTION_LENGTH)) %>%     # find the length of each section
  ungroup() %>%                                           # ungroup previous grouping by section number
  summarise(total_length = sum(section_length)) %>%       # sum section length of each section
  as.double()                                             # allows for 64 bit storage (increase precision with more significant digits)
    
total_spines     <- as.numeric(nrow(df))                  # number of rows = number of spines
density_overall  <- total_spines/total_length             # calculate density of dendritic segment
density_mushroom <- sum(as.numeric(df$TYPE == "mushroom"),# find total number of mushroom spines and
                        na.rm = TRUE)/total_length        # divide by total length to find mushroom density
density_thin     <- sum(as.numeric(df$TYPE == "thin"),    # as above for thin spines
                        na.rm = TRUE)/total_length        
density_stubby   <- sum(as.numeric(df$TYPE == "stubby"),  # as above for stubby spines
                        na.rm = TRUE)/total_length        
df               <- df[order(df$SOMA_DISTANCE), ]         # order data by SOMA_DISTANCE (attachment point)
data_coords      <- data.frame(df$X, df$Y, df$Z)          # create data frame with spine head coordinates


## Shane's original clustering analysis: uses soma distance (attachment point) and agnes
## agglomerative hierarchical clustering (UPGMA) for an 'ac' (agglomeration coefficient) to
## compare to a randomly generated population of ac's. This section of code also includes
## the original nearest neighbor analysis using euclidean distance of attachment points
    
agn <- agnes(df$SOMA_DISTANCE,                               # from package 'cluster'
              metric = "euclidean",                          # values are euclidean distances between soma distances
              method = "average")                            # "average" = UPGMA
    
dist_ac        <- as.matrix(dist(df$SOMA_DISTANCE))          # create distance matrix of soma distances, [1] is distance from itself or '0'
df$nn_dist_ac  <- apply(dist_ac, 2, function(x) sort(x)[2])  # find first nearest neighbor  (distance from 2nd closest) [2]
df$nn2_dist_ac <- apply(dist_ac, 2, function(x) sort(x)[3])  # find second nearest neighbor (distance from 3rd closest) [3]
df$nn3_dist_ac <- apply(dist_ac, 2, function(x) sort(x)[4])  # find third nearest neighboy  (distance from 4th closest) [4]
dendrite_ac    <- agn$ac                                     # define agglomeration coefficent for the dendrite
ac_test        <- list()                                     # initialize list for random ac's from next for-loop  
possible_dist  <- seq(1, total_length, by = 0.01)           # create list from 1 to total length of dendrite by 0.01 increments

for(k in 1:50) {  ## random ac for-loop, use random distances from soma to generate list of randomized ac values
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

## end Shane's original clustering analysis

## 1D clustering (from spine attachment points)

dist_1D          <- as.matrix(dist(df$SOMA_DISTANCE))                  # creates a distance matrix using distance from soma
UPGMA_obsv_1D    <- IdClusters(dist_1D,                                # runs cluster analysis on observed spines(cutoff from Yadav paper)
                               method = "UPGMA",                       # use UPGMA as clustering method
                               cutoff = 0.75,                          # cutoff from Yadav paper
                               showPlot = FALSE)                       # supress dendrogram
UPGMA_obsv_1D$rn <- rownames(UPGMA_obsv_1D)                            # add column with row names to keep spine ID (not sure if this step is necessary)
cluster_all_1D   <- cbind(UPGMA_obsv_1D, df$SOMA_DISTANCE)             # add cluster number, spine ID(row number) and corresponding distance from soma to same data frame
cluster_all_1D   <- cluster_all_1D[order(cluster_all_1D$cluster), ]    # orders data by cluster number
cluster_freq_1D  <- table(cluster_all_1D$cluster)                      # create table with # of spines in each cluster
cluster_freq_1D  <- as.data.frame(cluster_freq_1D)                     # convert above table to data frame
cluster_freq_1D$Freq         <- as.numeric(cluster_freq_1D$Freq)       # turn the "cluster" column to number
cluster_freq_1D$is_clustered <- as.numeric(cluster_freq_1D$Freq > 1)   # return 1 if  > 1 spine in cluster, 0 if not-- 1s reflect true "cluster" 
num_clusters_1D     <- sum(cluster_freq_1D$is_clustered, na.rm = TRUE) # count # of clusters on this segment
spines_clustered_1D <- sum(cluster_freq_1D$Freq > 1,     na.rm = TRUE) # sum total number of spines that are in clusters of > 1 spine
spines_clustered_1D <- as.matrix(spines_clustered_1D)                  # convert # of spines clustered to matrix
spines_clustered_1D <- as.numeric(spines_clustered_1D)                 # convert to numerical form
spines_clustered_1D[is.na(spines_clustered_1D)] <- 0                   # change any NA values to 0
spines_not_1D <- as.numeric(total_spines - spines_clustered_1D)        # total spines minus number of spines clustered = number of spines not clustered
test_spines_clustered_all_1D <- data.frame()                           # initialize dataframe for # of spines clustered for all random samples
test_spines_clustered_all_1D <- colnames("test spines clustered 1D")   # name the data frame column
test_num_clusters_all_1D     <- data.frame()                           # initialize dataframe for # of clusters for all random samples
test_num_clusters_all_1D     <- colnames("test number of clusters 1D") # name the data frame column

## 1D random spines for-loop:
for(j in 1:50){                                                        # generate list of randomly generated attachment points to calculate random clustering
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

  
  
# 3D clustering analysis  
  dist_3D <- as.matrix(dist.xyz(data_coords)) # creates distance matrix of spine head coordinates
  UPGMA_obsv_3D <- IdClusters(dist_3D, method = "UPGMA", cutoff = 0.75, showPlot = FALSE) #run cluster analysis with cutoff used in Yadav paper, gives cluster number associated with spine
  UPGMA_obsv_3D$rn <- rownames(UPGMA_obsv_3D)   #adds a column with rownows to keep spine ID, not sure if this step is neccessary
  cluster_all_3D <- cbind(UPGMA_obsv_3D, data_coords)   #get the coordinates for each spine ID, if kept in row name/number order, will correctly correspond to each spine
  cluster_all_3D <- cluster_all_3D[order(cluster_all_3D$cluster), ] #sort data frame by cluster number
  cluster_freq_3D <- table(cluster_all_3D$cluster) #create a table counting how many times each cluster Variable occurs (i.e. how many spines in each cluster)
  cluster_freq_3D <- as.data.frame(cluster_freq_3D)
  cluster_freq_3D$Freq <- as.numeric(cluster_freq_3D$Freq)
  cluster_freq_3D$is_clustered <- as.numeric(cluster_freq_3D$Freq > 1) # create column where 1 means there is more than one spine in a cluster or 0 if just 1
  num_clusters_3D <- sum(cluster_freq_3D$is_clustered, na.rm = TRUE) #count how many 1s to determine how many clusters (spines > 1) in the segment
  spines_clustered_3D <- cluster_freq_3D %>% 
    group_by(is_clustered) %>%
    filter(is_clustered == 0) %>% 
    summarise(spines_clustered_3D = sum(Freq, na.rm = TRUE))
  ifelse(is.null(spines_clustered_3D) ==TRUE, spines_clustered_3D <- 0, spines_clustered_3D <- spines_clustered_3D[1,2])
  spines_clustered_3D <- as.numeric(spines_clustered_3D) #converts to numerical form
  spines_not_3D <- as.numeric(total_spines - spines_clustered_3D)
  test_spines_clustered_all_3D <- data.frame()
  test_spines_clustered_all_3D <- colnames("test spines clustered 3D")
  test_num_clusters_all_3D <- data.frame()
  test_num_clusters_all_3D <- colnames("test number of clusters 3D")

  df$nn_dist_3D  <- apply(dist_3D, 2, function(x) sort(x)[2])  # find first nearest neighbor  (distance from 2nd closest) [2]
  df$nn2_dist_3D <- apply(dist_3D, 2, function(x) sort(x)[3])  # find second nearest neighbor (distance from 3rd closest) [3]
  df$nn3_dist_3D <- apply(dist_3D, 2, function(x) sort(x)[4])  # find third nearest neighboy  (distance from 4th closest) [4]
  
# 3D random spines for loop
  for(j in 1:50){
    test_data_X <- data.frame(sample(df$X), df$Y, df$Z) # randomize the X's, Y's, and Z's to make a "biologically plausible" dataframe.
    colnames(test_data_X) <- c( "x", "Y", "Z")
    test_data_Y <- data.frame(df$X, sample(df$Y), df$Z)
    colnames(test_data_Y) <- c( "x", "Y", "Z")
    test_data_Z <- data.frame(df$X, df$Y, sample(df$Z))
    colnames(test_data_Z) <- c( "x", "Y", "Z")
    test_data   <- rbind(test_data_X, test_data_Y, test_data_Z)
    test_data_final <- data.frame(sample_n(test_data, total_spines))
    test_dist_3D <- as.matrix(dist(test_data_final)) #creates distance matrix for random sample

    UPGMA_test_3D <- IdClusters(test_dist_3D, method = "UPGMA", cutoff = 0.75, showPlot = FALSE) #run cluster analysis with cutoff used in Yadav paper
    UPGMA_test_3D$rn <- rownames(UPGMA_test_3D)   #adds a column with rownows to keep spine ID, not sure if this step is neccessary
    cluster_all_test_3D <- cbind(UPGMA_test_3D, test_data_final)   #get the coordinates for each spine ID, if kept in row name/number order, will correctly correspond to each spine
    cluster_all_test_3D <- cluster_all_test_3D[order(cluster_all_test_3D$cluster), ] #sort data frame by cluster number
    cluster_freq_test_3D <- table(cluster_all_test_3D$cluster) #create a table counting how many times each cluster Variable occurs (i.e. how many spines in each cluster)
    cluster_freq_test_3D <- as.data.frame(cluster_freq_test_3D)
    cluster_freq_test_3D$is_clustered <- as.numeric(cluster_freq_test_3D$Freq > 1) # create column where 1 means there is more than one spine in a cluster or 0 if just 1
    num_clusters_test_3D <- sum(cluster_freq_test_3D$is_clustered, na.rm = TRUE)
    num_clusters_test_3D[is.na(num_clusters_test_3D)] <- 0
    num_clusters_test_3D <- data.frame(num_clusters_test_3D)
    test_num_clusters_all_3D <- rbind(test_num_clusters_all_3D, num_clusters_test_3D)
    spines_clustered_test_3D <- cluster_freq_test_3D %>% 
      group_by(is_clustered) %>% 
      filter(is_clustered == 1) %>% 
      summarise(spines_clustered_test_3D = sum(Freq))
    ifelse(is.null(spines_clustered_test_3D) == TRUE, spines_clustered_test_3D <- 0, spines_clustered_test_3D <- spines_clustered_test_3D[1,2])
    spines_not_test_3D <- as.numeric(total_spines - spines_clustered_test_3D)
    spines_clustered_test_3D[is.na(spines_clustered_test_3D)] <- 0
    spines_clustered_test_3D <- data.frame(spines_clustered_test_3D)
    test_spines_clustered_all_3D <- rbind(test_spines_clustered_all_3D, spines_clustered_test_3D)
    
  } # end of random spines 3D for-loop
  
  
  test_spines_clustered_all_3D <- as.numeric(as.matrix(test_spines_clustered_all_3D))
  std_test_3D <- sd(test_spines_clustered_all_3D)
  mean_test_3D <- mean(test_spines_clustered_all_3D)
  curve_dnorm_3D <- dnorm(test_spines_clustered_all_3D, mean_test_3D, std_test_3D)
  std_curve_3D <- sd(curve_dnorm_3D)
  mean_curve_3D <- mean(curve_dnorm_3D)
  Cscore_3D <- pnorm(spines_clustered_3D, mean_test_3D, std_test_3D)
  
  
#adding data to running list in data_all file
  df$dendrite_ac         <- dendrite_ac
  df$Cscore_1D           <- Cscore_1D
  df$Cscore_3D           <- Cscore_3D
  df$num_clusters_3D     <- num_clusters_3D
  df$spines_clustered_3D <- spines_clustered_3D
  df$spines_not_3D       <- spines_not_3D
  df$density_overall     <- density_overall # add these data to df
  df$density_mushroom    <- density_mushroom
  df$density_thin        <- density_thin
  df$density_stubby      <- density_stubby
  data_all               <- rbind(data_all, df) # adds df as next row in the data_all file
  
} # end of file for-loop 




#For following comments 1-6 we are applying the same code to two separate data-frames of different lengths
#data_all contains all values for each spine (row), which means that some values that would be the same for a group of spines in the same file
#(like overall density, or overall Cscore), are repeated for each row
#data_all_by_file groups these values by file name so there are no repeats in values when averaging density, Cscore, etc.

data_all_by_file <- data_all %>% 
  group_by(file) %>% 
  summarise(Cscore_1D           = mean(Cscore_1D),
            cScore_ac_1D        = mean(cScore_ac_1D),
            dendrite_ac         = mean(dendrite_ac),
            num_clusters_3D     = mean(num_clusters_3D), 
            spines_clustered_3D = mean(spines_clustered_3D), 
            spines_not_3D       = mean(spines_not_3D), 
            Cscore_3D           = mean(Cscore_3D), 
            density_overall     = mean(density_overall), 
            density_mushroom    = mean(density_mushroom), 
            density_thin        = mean(density_thin), 
            density_stubby      = mean(density_stubby),
            dendrite_ac         = mean(dendrite_ac))


#1) Add column with group to master data file
#2) Pull animal number ID out of file name, creates column with animal number, also removes the L/N following animal number in file name
#3) Pulls off the "L" or "N" from file ID and add a column, changes to "labeled" or "not labeled", respectively

data_all$group         <- group_name
data_all$animal_num <- lapply(data_all$file,
                              function(x) unlist(strsplit(x, "-"))[2]) #pulls animal ID out of file name, with retro label attached

data_all$retro_label         <- substring(data_all$animal_num, # adds column reported retro-label
                                          nchar(data_all$animal_num),
                                          nchar(data_all$animal_num))

data_all$animal_num <- lapply(data_all$animal_num, 
                              function(x) unlist(strsplit(x, "L"))[1]) #removes L or N from animal number

data_all$animal_num <- lapply(data_all$animal_num, 
                              function(x) unlist(strsplit(x, "N"))[1]) #removes L or N from animal number

data_all$retro_label <- replace(data_all$retro_label, 
                                data_all$retro_label == "L", "labeled") #adds labeled for L and not labeled for N
data_all$retro_label <- replace(data_all$retro_label, 
                                data_all$retro_label == "N", "not labeled")

data_all_by_file$group <- group_name
data_all_by_file$animal_num <- lapply(data_all_by_file$file,           #pulls animal ID out of file name, with retro label attached
                                      function(x) unlist(strsplit(x, "-"))[2])
data_all_by_file$retro_label <- substring(data_all_by_file$animal_num, # adds column reported retro-label
                                          nchar(data_all_by_file$animal_num), 
                                          nchar(data_all_by_file$animal_num))
data_all_by_file$animal_num <- lapply(data_all_by_file$animal_num, 
                                      function(x) unlist(strsplit(x, "L"))[1])#removes L or N from animal number
data_all_by_file$animal_num <- lapply(data_all_by_file$animal_num, 
                                      function(x) unlist(strsplit(x, "N"))[1])#removes L or N from animal number

data_all_by_file$retro_label <- replace(data_all_by_file$retro_label, 
                                        data_all_by_file$retro_label == "L", #adds labeled for L and not labeled for N
                                        "labeled") 
data_all_by_file$retro_label <- replace(data_all_by_file$retro_label, 
                                        data_all_by_file$retro_label == "N", 
                                        "not labeled")


#4) Pulls of dendrite location from file name and adds to column, changes letters to full location names(basal, proximal, distal)
data_all$PDB <- lapply(data_all$file, function(x) unlist(strsplit(x, "-"))[3])
data_all$PDB <- replace(data_all$PDB, data_all$PDB == "p", "prox")
data_all$PDB <- replace(data_all$PDB, data_all$PDB == "d", "dist")
data_all$PDB <- replace(data_all$PDB, data_all$PDB == "b", "basal")

data_all_by_file$PDB <- lapply(data_all_by_file$file, 
                               function(x) unlist(strsplit(x, "-"))[3]) # pulls dendrite location out of name and adds column
data_all_by_file$PDB <- replace(data_all_by_file$PDB, 
                                data_all_by_file$PDB == "b", "basal")
data_all_by_file$PDB <- replace(data_all_by_file$PDB, 
                                data_all_by_file$PDB == "p", "prox")
data_all_by_file$PDB <- replace(data_all_by_file$PDB, 
                                data_all_by_file$PDB == "d", "dist")


#5) Pulls day of imaging from file name (letter at the end) and add a column
data_all$stack <- lapply(data_all$file, 
                         function(x) unlist(strsplit(x, "-"))[4])
data_all$stack <- unlist(data_all$stack)

data_all_by_file$stack <- lapply(data_all_by_file$file, 
                                 function(x) unlist(strsplit(x, "-"))[4]) 
data_all_by_file$stack <- unlist(data_all_by_file$stack)


#6) Change animal_num and dendrite location(PDB) to a vector rater than a list, which is important for executing upcoming code
data_all$animal_num <- unlist(data_all$animal_num)
data_all$PDB        <- unlist(data_all$PDB)

data_all_by_file$animal_num <- unlist(data_all_by_file$animal_num)
data_all_by_file$PDB        <- unlist(data_all_by_file$PDB)

# Change columns to numeric
data_all$RAYBURST_VOLUME <- as.numeric(data_all$RAYBURST_VOLUME)
data_all$MAX_DTS         <- as.numeric(data_all$MAX_DTS)



#Summary Data from data_all_by_file
location_avg_file <- data_all_by_file %>% 
  group_by(group, animal_num, retro_label, PDB) %>% 
  summarise(num_clusters_3D     = mean(num_clusters_3D),
            spines_clustered_3D = mean(spines_clustered_3D),
            spines_not_3D       = mean(spines_not_3D),
            Cscore_3D           = mean(Cscore_3D),
            density_overall     = mean(density_overall),
            density_mushroom    = mean(density_mushroom),
            density_thin        = mean(density_thin),
            density_stubby      = mean(density_stubby))

animal_labeled_avg_file <- data_all_by_file %>% 
  group_by(group, animal_num, retro_label) %>% 
  summarise(num_clusters_3D     = mean(num_clusters_3D),
            spines_clustered_3D = mean(spines_clustered_3D),
            spines_not_3D       = mean(spines_not_3D),
            Cscore_3D           = mean(Cscore_3D),
            density_overall     = mean(density_overall),
            density_mushroom    = mean(density_mushroom),
            density_thin        = mean(density_thin),
            density_stubby      = mean(density_stubby))

#Summary data for data_all
spine_type_ave <- data_all %>% 
  group_by(group, animal_num, retro_label, TYPE) %>% 
  summarise(spine_vol_overall    = mean(RAYBURST_VOLUME, na.rm = TRUE), 
            spine_length_overall = mean(MAX_DTS, na.rm = TRUE))

spine_type_by_location_ave <- data_all %>% 
  group_by(group, animal_num, stack, retro_label, PDB, TYPE) %>% 
  summarise(spine_vol_overall    = mean(RAYBURST_VOLUME, na.rm = TRUE), 
            spine_length_overall = mean(MAX_DTS, na.rm = TRUE))

location_avg <- data_all %>% 
  group_by(group, animal_num, retro_label, PDB) %>% 
  summarise(spine_vol_overall    = mean(RAYBURST_VOLUME, na.rm = TRUE), 
            spine_length_overall = mean(MAX_DTS, na.rm = TRUE))

animal_by_label_ave <- data_all %>% 
  group_by(group, animal_num, retro_label) %>% 
  summarise(spine_vol_overall    = mean(RAYBURST_VOLUME, na.rm = TRUE), 
            spine_length_overall = mean(MAX_DTS, na.rm = TRUE))


nn_dist_ac <- data_all$nn_dist_ac
nn_dist_ac <- cbind(nn_dist_ac, data_all$nn2_dist_ac)
nn_dist_ac <- cbind(nn_dist_ac, data_all$nn3_dist_ac)
nn_dist_ac_colnames <- c("nn1_1D", "nn2_1D", "nn3_1D")
colnames(nn_dist_ac) <- nn_dist_ac_colnames

nn_dist_3D <- data_all$nn_dist_3D
nn_dist_3D <- cbind(nn_dist_3D, data_all$nn2_dist_3D)
nn_dist_3D <- cbind(nn_dist_3D, data_all$nn3_dist_3D)
nn_dist_3D_colnames <- c("nn1_3D", "nn2_3D", "nn3_3D")
colnames(nn_dist_3D) <- nn_dist_3D_colnames


dir.create("analysis5")  #creates a directory to create a file on the computer
savePath <- paste(filePath, "/", "analysis5", collapse = "/", sep = "") # creates the path where files can be saved
setwd(savePath) # sets working directory to the path created above
write.csv(data_all, "data_all.csv")
write.csv(data_all_by_file, "data_all_by_file.csv")
write.csv(location_avg_file, "location_avg_file.csv")
write.csv(animal_labeled_avg_file, "animal_labeled_avg_file.csv")
write.csv(spine_type_ave, "spine_type_ave.csv")
write.csv(spine_type_by_location_ave, "spine_type_by_location_ave.csv")
write.csv(animal_by_label_ave, "animal_by_label_ave.csv")
write.csv(nn_dist_ac, "nn_1D.csv")
write.csv(nn_dist_3D, "nn_3D.csv")

time_elapsed_hours <- (timer/1000)/60/60
print(time_elapsed_hours)
winDialog("ok", "Code complete")

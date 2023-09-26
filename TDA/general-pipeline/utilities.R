
#Some functions here are from Peter Bubenik's R-labs that can be found at 
#https://people.clas.ufl.edu/peterbubenik/intro-to-tda/

get_csvs <- function(path) {
    # get list of discretized CSVs recursively
    files <- list.files(path, pattern="\\.csv$", include.dirs=TRUE, recursive=TRUE)
}


concat_path <- function(path, filename) {
  # make sure path has file separator at the end
  n <- nchar(path)
  if (substr(path, n, n) != .Platform$file.sep) {
    path <- paste0(path, .Platform$file.sep)
  }

  # return correct joined path to file
  return(paste0(path, filename))
}


# remove dimension zero points from a persistence diagram
remove_dimzero_from_diagrams <- function(PD){
  zeroes <- which(PD[, 1] == 0)
  new_PD <- PD[-zeroes,] 
  return(new_PD)
}


#compute and save a list of persistence diagrams
compute_PDs <- function(csv_files_path, group_name, save_location, types, index = NaN){
  
  # recursively find CSV files within directory and return file-names in a list
  data_files <- get_csvs(csv_files_path)
  
  PDs <- list()
  max_birth <- list() 
  max_death <- list() 
  for (i in 1:length(data_files)){
    print(sprintf("Processing file %s", data_files[i]))
    file_path <- concat_path(csv_files_path, data_files[i])
    cells <- read.csv(file_path)
    # if no index provided, use last column
    if (is.na(index)) index <- ncol(cells)
    
    # if there are cell types specified, only include those rows
    if (!any(is.na(types))) {
      if (length(types) > 1) {
        cells <- cells[cells[,index] %in% types,]    # if a vector is provided
      } else {
        cells <- cells[cells[,index] == types,]    # if a single value
      }
    }
    
    #compute persistence homology using Delaunay complex filtration (also known as Alpha complex filtration)
    PH <-  alphaComplexDiag(cells[,1:2], maxdimension = 1, library = c("GUDHI", "Dionysus"), location = TRUE)
    PDs[[i]] <- PH[["diagram"]]
    PDs[[i]] <- remove_dimzero_from_diagrams(PDs[[i]])
    PDs[[i]][,2:3] <- sqrt(PDs[[i]][,2:3])
    #max birth of 1-degree persistence diagrams 
    max_birth[[i]] <- max(PDs[[i]][,2])
    #max death of 1-degree persistence diagrams 
    max_death[[i]] <- max(PDs[[i]][,3])
  }
  
  max_birth <- max(unlist(max_birth))
  max_death <- max(unlist(max_death))
  
  #save max birth and max death to max-values.csv file
  df <- read.csv(concat_path(save_location,"max-values.csv"))
  max_values <- rbind(df, data.frame(maxbirth = max_birth, maxdeath = max_death))
  write.csv(max_values, concat_path(save_location,"max-values.csv"), row.names=FALSE)
  
  path <- concat_path(save_location, group_name)
  
  #save persistence diagrams
  save(PDs, file=concat_path(path, paste0(group_name,"_PDs.RData")))
}


#plot persistence diagrams from a list
plot_PDs <- function(group_name, birth, death, save_location, save_plot){
  
  #location of saved PDs
  path <- concat_path(save_location, group_name)
  PD_list <- get(load(concat_path(path, paste0(group_name, "_PDs.RData"))))
  
  par(pty="s")
  max_radius <- max(birth, death)+5
  for (i in 1:length(PD_list)){
      plot(PD_list[[i]][,2:3], asp=1, xlim=c(0, max_radius) , ylim=c(0, max_radius), xlab='', ylab='', col='orange', bty="o", pch=20, cex=1)
      abline(0,1)
      
      #if TRUE save plot
      if (save_plot){
        dev.copy(pdf, concat_path(path, paste0(group_name,"_PD_", i, ".pdf")))
        invisible(dev.off())
      }
  }
}


#plot representative cycles that persist (live) over a certain threshold
plot_representative_cycles <- function(csv_files_path, cell_types, threshold, group_name, save_location, save_plot){
  
  # recursively find CSV files within directory and return file-names in a list
  data_files <- get_csvs(csv_dir_path)
  
  par(pty="s")
  for (i in 1:length(data_files)){
    print(sprintf("Processing file %s", data_files[i]))
    file_path <- concat_path(csv_files_path, data_files[i])
    cells <- read.csv(file_path)
    cells <- cells[which(is.element(cells[,8], cell_types)),]
    #compute persistence homology using Delaunay complex filtration (also known as Alpha complex filtration)
    filtration <- alphaComplexFiltration(cells[,1:2], printProgress = TRUE)
    PH <-  alphaComplexDiag(cells[,1:2], maxdimension = 1, library = c("GUDHI", "Dionysus"), location = TRUE)
    PD <- PH[["diagram"]]
    #plot cycles that persist over specific persistence_param
    ones <- which(PD[,1] == 1)
    if (length(ones) > 1 ){
      plot(filtration[["coordinates"]], pch = 19, cex=0.1, ylab="", xlab="", axes = FALSE)
      for (m in ones[1]:(length(ones)+ones[1]-1)){
        cycles <- PH[["cycleLocation"]][m]
        if ((sqrt(PD[m,3]) - sqrt(PD[m,2])  > threshold) ){
          for (s in 1:length(cycles)){
            for (l in 1:dim(cycles[[s]])[1]){
              lines(cycles[[s]][l,,], col="blue", lwd=1)
            }
          } 
        } 
      }
    }
    #location of saved PDs
    path <- concat_path(save_location, group_name)
    
    #if TRUE save plot
    if (save_plot){
      dev.copy(pdf, concat_path(path, paste0(group_name,"_repr_cycles_", i, ".pdf")))
      invisible(dev.off())
    }
  }
}


#compute and save a list of persistence landscapes
compute_PLs <- function(group_name, birth, death, discr_step, save_location){
  
  #location of saved PLs
  path <- concat_path(save_location, group_name)
  PD_list <- get(load(concat_path(path, paste0(group_name, "_PDs.RData"))))
  
  max_x <- (death+birth)/2 + 5
  radius_values <- seq(0, max_x, discr_step)
  PLs <- list()
  for (i in 1:length(PD_list)){ 
    PLs[[i]] <- t(landscape(PD_list[[i]],dimension=1,KK=1:1000,radius_values))
  }
  
  #save persistence landscapes
  save(PLs, file=concat_path(path, paste0(group_name,"_PLs.RData")))
}


#plot a list of persistence landscapes
plot_PLs <- function(group_name, birth, death, discr_step, save_location, save_plot){
  
  #location of saved PLs
  path <- concat_path(save_location, group_name)
  PL_list <- get(load(concat_path(path, paste0(group_name, "_PLs.RData"))))
  
  par(pty="s")
  mycolors <- c("midnightblue","slateblue", "royalblue2", "deepskyblue4", "deepskyblue3",
                "turquoise4", "chartreuse4","olivedrab4", "olivedrab3", "yellowgreen", 
                "yellow3","gold3", "gold2","gold1", "gold" )
  
  max_x <- (death+birth)/2 + 5
  max_height <- death/2 
  radius_values <- seq(0, max_x, discr_step)
  for (i in 1:length(PL_list)){ 
    plot(radius_values, PL_list[[i]][1,], type="l", ylab="", xlab="", xlim=c(0, max_x) , ylim=c(0, max_height), ann=FALSE, bty="o",col=mycolors[1])
    for(k in 2:dim(PL_list[[i]])[1]){
      lines(radius_values,PL_list[[i]][k,],type="l",col=mycolors[k %% 15])
    }
    
    #if TRUE save plot
    if (save_plot){
      dev.copy(pdf, concat_path(path, paste0(group_name,"_PL_", i, ".pdf")))
      invisible(dev.off())
    }
  }
}


# convert a vector to a persistence landscape
landscape_from_vector <- function(PL_vector, radius_values){
  
  m <- length(radius_values)
  K <- length(PL_vector)/m
  PL <- Matrix(0, nrow = K, ncol=m, sparse = TRUE)
  for (i in 1:K){
    PL[i,1:m] <- PL_vector[(1+(i-1)*m):(i*m)]
  }
  return(PL)
}


# Matrix of persistence landscape row vectors from list of persistence landscapes
landscape_matrix_from_list <- function(PL_list){
  
  n <- length(PL_list)
  m <- ncol(PL_list[[1]])
  max_depth <- integer(n)
  for (i in 1:n)
    max_depth[i] <- nrow(PL_list[[i]])
  K <- max(max_depth)
  PL_matrix <- Matrix(0, nrow = n, ncol = m*K, sparse = TRUE)
  for (i in 1:n)
    for (j in 1:max_depth[i])
      PL_matrix[i,(1+(j-1)*m):(j*m)] <- PL_list[[i]][j,]
  return(PL_matrix)
}


#compute and save the average persistence landscape
compute_avgPL <- function(group_name, birth, death, discr_step, save_location){
  
  #location of saved PLs
  path <- concat_path(save_location, group_name)
  PL_list <- get(load(concat_path(path, paste0(group_name, "_PLs.RData"))))
  
  max_x <- (death+birth)/2 + 5
  radius_values <- seq(0, max_x, discr_step)
  PL_matrix <- landscape_matrix_from_list(PL_list)
  average_PL_vector <- colMeans(PL_matrix, sparseResult = TRUE)
  average_PL <- landscape_from_vector(average_PL_vector, radius_values)
  
  #save the list of persistence diagrams, the max birth and max death values
  save(average_PL, file=concat_path(path, paste0(group_name,"_avgPL.RData")))
}


#plot average persistence landscape
plot_avgPL <- function(group_name, birth, death, discr_step, save_location, save_plot){
  
  #location of a saved average PL
  path <- concat_path(save_location, group_name)
  average_PL <- get(load(concat_path(path, paste0(group_name, "_avgPL.RData"))))
  
  par(pty="s")
  mycolors <- c("midnightblue","slateblue", "royalblue2", "deepskyblue4", "deepskyblue3",
                "turquoise4", "chartreuse4","olivedrab4", "olivedrab3", "yellowgreen", 
                "yellow3","gold3", "gold2","gold1", "gold" )
  
  max_x <- (death+birth)/2 + 5
  radius_values <- seq(0, max_x, discr_step)
  max_height <- death/2 
  plot(radius_values, average_PL[1,], type="l", ylab="", xlab="", xlim=c(0, max_x) , ylim=c(0, max_height), ann=FALSE, bty="o",col=mycolors[1])
  for(k in 2:dim(average_PL)[1]){
    lines(radius_values,average_PL[k,],type="l",col=mycolors[k %% 15])
  }
  
  #if TRUE save plot
  if (save_plot){
    dev.copy(pdf, concat_path(path, paste0(group_name, "_avgPL_plot.pdf")))
    invisible(dev.off())
  }
}


#compute the difference between two average persistence landscapes 
plot_avgPLs_difference <- function(groups, birth, death, discr_step, save_location, save_plot){
  
  #location of the saved first average PL
  path_1 <- concat_path(save_location, groups[1])
  average_PL1 <- get(load(concat_path(path_1, paste0(groups[1], "_avgPL.RData"))))
  
  #location of the saved second average PL
  path_2 <- concat_path(save_location, groups[2])
  average_PL2 <- get(load(concat_path(path_2, paste0(groups[2], "_avgPL.RData"))))
  
  max_x <- (max_birth+max_death)/2 + 5
  max_height <- max_death/2
  radius_values <- seq(0, max_x, discr_step)
  
  difference_matrix <- average_PL1 - average_PL2
  
  par(pty="s")
  mycolors <- c("midnightblue","slateblue", "royalblue2", "deepskyblue4", "deepskyblue3",
                "turquoise4", "chartreuse4","olivedrab4", "olivedrab3", "yellowgreen", 
                "yellow3","gold3", "gold2","gold1", "gold" )
  
  plot(radius_values, difference_matrix[2,], type="l", ylab="", xlab="", xlim=c(0, max_x), ylim=c(-max_height, max_height), ann=FALSE, bty="o",col=mycolors[1])
  for(k in 2:dim(difference_matrix)[1]){
    lines(radius_values, difference_matrix[k,],type="l",col=mycolors[k %% 15])
  }
  
  #if TRUE save plot
  if (save_plot){
    dev.copy(pdf, concat_path(save_location, "avgPL_difference.pdf"))
    invisible(dev.off())
  }
}


#Euclidean distance between two vectors 
euclidean_distance <- function(u, v) sqrt(sum((u - v) ^ 2))


#permutation test for PLs
permutation_test_for_PLs <- function(groups, save_location, nrepeats){
  
  #location of saved PLs from group1
  path_1 <- concat_path(save_location, groups[1])
  PL1 <- get(load(concat_path(path_1, paste0(groups[1], "_PLs.RData"))))
  
  #location of saved PLs from group2
  path_2 <- concat_path(save_location, groups[2])
  PL2 <- get(load(concat_path(path_2, paste0(groups[2], "_PLs.RData"))))
  
  PL1_matrix <- landscape_matrix_from_list(PL1)
  PL2_matrix <- landscape_matrix_from_list(PL2)
  # append zeros if necessary so that the matrices have the same number of columns
  num_columns <- max(ncol(PL1_matrix),ncol(PL2_matrix))
  PL1_matrix <- cbind(PL1_matrix, Matrix(0,nrow=nrow(PL1_matrix),ncol=num_columns-ncol(PL1_matrix)))
  PL2_matrix <- cbind(PL2_matrix, Matrix(0,nrow=nrow(PL2_matrix),ncol=num_columns-ncol(PL2_matrix)))
  
  M <- rbind(PL1_matrix, PL2_matrix)
  n <- nrow(M) # number of total data points
  k <- nrow(PL1_matrix) # number of data points in first group
  observation <- euclidean_distance(colMeans(PL1_matrix),colMeans(PL2_matrix))
  ngreater_distances <- 0
  for (i in 1:nrepeats){
    permutation <- sample(1:n)
    distance <- euclidean_distance(colMeans(M[permutation[1:k],]),colMeans(M[permutation[(k+1):n],]))
    if (distance >= observation)
      ngreater_distances <- ngreater_distances + 1
  }
  
  #compute p-value
  p_value <- ngreater_distances/nrepeats
  if (p_value == 0){
    p_value <- 1/nrepeats #include a permutation that gives an observed value
  }
  
  return(print(sprintf("p-value %.4f:", p_value)))
}


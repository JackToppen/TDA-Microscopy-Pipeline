library(ggplot2)
library(TDA)


make_PDs <- function(csv_path, out_path, types, index=NA) {
  # hold max birth/death and features for plotting and landscapes
  max_birth <- 0
  max_death <- 0
  max_features <- 0
  
  # iterate through list of CSVs with locations to generate diagrams
  files <- get_files(csv_path, "\\.csv$")
  for (file in files) {
    # read the CSV to a data frame
    csv_file <- concat_paths(csv_path, file)
    df <- read.csv(csv_file)
    
    # if no index provided, use last column
    if (is.na(index)) index <- ncol(df)
  
    # if there are cell types specified, only include those rows
    if (!any(is.na(types))) {
      if (length(types) > 1) {
        df <- df[df[,index] %in% types,]    # if a vector is provided
      } else {
        df <- df[df[,index] == types,]    # if a single value
      }
    }
    
    # compute persistence homology using Delaunay complex filtration (aka Alpha complex filtration)
    out <- alphaComplexDiag(df[,1:2], maxdimension=1, library=c("GUDHI", "Dionysus"), location=TRUE)
    PD <- out[["diagram"]]
    
    # remove dimension zero points
    PD <- PD[PD[,1] != 0,,drop=F]
    
    # only continue if features exist
    if (nrow(PD) > 0) {
      PD[,2:3] <- sqrt(PD[,2:3])
      
      # hold max birth/death all 1-degree persistence diagrams
      if (nrow(PD) > 1) {
        birth <- max(PD[,2])
        death <- max(PD[,3])
      } else {
        birth <- PD[1,2]
        death <- PD[1,3]
      }
      if (birth > max_birth) max_birth <- birth
      if (death > max_death) max_death <- death
      
      # additionally get max number of features
      features <- nrow(PD)
      if (features > max_features) max_features <- features
      
      # get path to diagram CSV output
      no_ext <- sub("\\.csv$", "", file)
      filename <- paste0(no_ext, .Platform$file.sep, basename(no_ext), "_PD.csv")
      PD_path <- concat_paths(out_path, filename)
      
      # make directory to output and write the CSV
      make_dir(dirname(PD_path))
      write.csv(PD, PD_path, row.names=FALSE)
    }
  }
   
  # save max birth/death for plotting and landscape generating
  rds_file <- concat_paths(out_path, "maxes.rds")
  saveRDS(c(max_birth, max_death, max_features), file=rds_file)
}


plot_PDs <- function(out_path, color, resolution, plot_max=NA) {
  # read the max/birth RDS file
  rds_file <- concat_paths(out_path, "maxes.rds")
  maxes <- readRDS(rds_file)
  
  # iterate through all PD files
  files <- get_files(out_path, "\\PD.csv$")
  for (file in files) {
    # read the CSV and get diagram
    PD_file <- concat_paths(out_path, file)
    PD <- read.csv(PD_file)
    
    # convert filename from CSV to PNG and get path to new PNG
    local_path <- sub("\\.csv$", ".png", file)
    image_path <- concat_paths(out_path, local_path)
    
    # plot the persistence diagram (include extra space beyond max death)
    if (!is.na(plot_max)) {
      plot_lim <- plot_max
    } else {
      plot_lim <- 1.1 * maxes[2]
    }
    p <- ggplot(PD, aes(x=Birth, y=Death)) +
      geom_point(shape=16, color=color) +
      scale_x_continuous(expand=c(0, 0), limits=c(0, plot_lim)) + 
      scale_y_continuous(expand=c(0, 0), limits=c(0, plot_lim)) + 
      theme_classic() + xlab("birth") + ylab("death")
    
    # add a line with slope=1
    p + geom_abline(intercept=0, slope=1)
    
    # convert width to inches and write image file
    width <- resolution / 300    # 300 pixels per inch default
    ggsave(image_path, width=width, height=width)
  }
}


make_PLs <- function(out_path, percent=NA, n=NA, disc_step=NA) {
  # read the max/birth RDS file
  rds_file <- concat_paths(out_path, "maxes.rds")
  maxes <- readRDS(rds_file)
  
  # if no user-specified discretization step, base the step on the max death
  if (is.na(disc_step)) disc_step = maxes[2] / 200
  
  # PL limit for plotting/comparing
  x_max <- (maxes[1]+maxes[2])/2
  values <- seq(0, x_max, disc_step)
  
  # determine the number of landscape functions to calculate
  if (is.na(n)) {
    if (is.na(percent)) {
      n <- 1000    # default value
    } else {
      # otherwise compute percentage of max possible number of landscapes
      n <- as.integer(percent * maxes[3])    
    }
  }
  
  # iterate through all PD files
  files <- get_files(out_path, "\\PD.csv$")
  for (file in files) {
    # read the CSV and get diagram
    PD_file <- concat_paths(out_path, file)
    PD <- read.csv(PD_file)
    
    # generate the landscape (each column is a landscape function)
    PL <- as.data.frame(landscape(PD, dimension=1, KK=1:n, values))
    colnames(PL) <- paste0("func_", 1:n)    # name columns
    PL <- cbind(values, PL)    # add x-axis values
    
    # get path to landscape CSV output
    local_path <- sub("\\_PD.csv$", "_PL.csv", file)
    PL_file <- concat_paths(out_path, local_path)
    
    # write the CSV
    write.csv(PL, PL_file, row.names=FALSE)
  }
}
    

plot_PLs <- function(out_path, n, resolution) {
  # read the max/birth RDS file
  rds_file <- concat_paths(out_path, "maxes.rds")
  maxes <- readRDS(rds_file)
  
  # get plot limits
  x_max <- (maxes[1]+maxes[2])/2
  y_max <- maxes[2]/2
  
  # iterate through all PL files
  files <- get_files(out_path, "\\PL.csv$")
  for (file in files) {
    # read the CSV and get landscape
    PL_file <- concat_paths(out_path, file)
    PL <- read.csv(PL_file)
    
    # convert filename from CSV to PNG and get path to new PNG
    local_path <- sub("\\.csv$", ".png", file)
    image_path <- concat_paths(out_path, local_path)
    
    # generate and save persistence landscape plot
    plot_PL(PL, image_path, c(0, x_max), c(0, y_max), n, resolution)
  }
}


make_avgPL <- function(path, groups) {
  # generate average landscape for each group
  for (dir in groups) {
    # get all PL files in the specified directory
    dir_path <- concat_paths(path, dir)
    files <- get_files(dir_path, "\\PL.csv$")
    n <- length(files)
    
    # read the first PL to start average
    PL_file <- concat_paths(dir_path, files[[1]])
    total <- data.matrix(read.csv(PL_file))
    
    # iterate through all files
    for (i in 2:n) {
      # read the CSV and add landscape to total
      PL_file <- concat_paths(dir_path, files[[i]])
      total <- total + data.matrix(read.csv(PL_file))
    }
    
    # get path to landscape CSV output
    local_path <- paste0(basename(dir), "_avgPL.csv")
    avgPL_file <- concat_paths(dir_path, local_path)

    # get average PL and write the CSV
    avgPL <- total/n
    write.csv(avgPL, avgPL_file, row.names=FALSE)
  }
}


plot_avgPL <- function(path, groups, n, resolution) {
  # read the max/birth RDS file
  rds_file <- concat_paths(path, "maxes.rds")
  maxes <- readRDS(rds_file)
  
  # get plot limits
  x_max <- (maxes[1]+maxes[2])/2
  y_max <- maxes[2]/2
  
  # plot average landscape for each group
  for (dir in groups) {
    # get path to average landscape file
    dir_path <- concat_paths(path, dir)
    local_path <- paste0(basename(dir), "_avgPL.csv")
    
    # read the CSV and get landscape
    avgPL_file <- concat_paths(dir_path, local_path)
    avgPL <- read.csv(avgPL_file)
    
    # convert filename and save persistence landscape plot
    image_path <- sub("\\.csv$", ".png", avgPL_file)
    plot_PL(avgPL, image_path, c(0, x_max), c(0, y_max), n, resolution)
  }
}


plot_avgPL_diff <- function(path, groups, n, resolution) {
  # read the max/birth RDS file
  rds_file <- concat_paths(path, "maxes.rds")
  maxes <- readRDS(rds_file)
  
  # get plot limits
  x_max <- (maxes[1]+maxes[2])/2
  y_max <- maxes[2]/2
  
  # get all possible pairs and iterate through them
  pairs <- combn(groups, 2)
  for (i in 1:ncol(pairs)) {
    dir_1 <- pairs[1,i]
    dir_2 <- pairs[2,i]
    
    # get path to first average landscape file
    dir_path_1 <- concat_paths(path, dir_1)
    local_path_1 <- paste0(basename(dir_1), "_avgPL.csv")
    avgPL_file_1 <- concat_paths(dir_path_1, local_path_1)
    
    # get path to second average landscape file
    dir_path_2 <- concat_paths(path, dir_2)
    local_path_2 <- paste0(basename(dir_2), "_avgPL.csv")
    avgPL_file_2 <- concat_paths(dir_path_2, local_path_2)
    
    # read average persistence landscapes
    avgPL_1 <- data.matrix(read.csv(avgPL_file_1))
    avgPL_2 <- data.matrix(read.csv(avgPL_file_2))
    diffPL <- avgPL_1 - avgPL_2
    diffPL[,1] <- avgPL_1[,1]
    
    # change the file separator to dash if path contains it
    rename_1 <- sub(.Platform$file.sep, "-", dir_1)
    rename_2 <- sub(.Platform$file.sep, "-", dir_2)
    
    # plot difference of landscapes
    image_path <- concat_paths(path, paste0(rename_1, "_", rename_2, "_diff.png"))
    plot_PL(diffPL, image_path, c(0, x_max), c(-y_max, y_max), n, resolution)
  }
}


perm_test <- function(path, groups, depth=30, reps=10000, all=FALSE) {
  # get all possible pairs and iterate through them
  pairs <- combn(groups, 2)
  for (i in 1:ncol(pairs)) {
    dir_1 <- pairs[1,i]
    dir_2 <- pairs[2,i]
    
    # get a matrix with rows as each PL within the first group
    dir_path_1 <- concat_paths(path, dir_1)
    PL_mat_1 <- PLs_to_mat(dir_path_1, depth, all)
    
    # get a matrix with rows as each PL within the second group
    dir_path_2 <- concat_paths(path, dir_2)
    PL_mat_2 <- PLs_to_mat(dir_path_2, depth, all)
    
    # compute average PL vector from each matrix
    avgPL_1 <- colMeans(PL_mat_1)
    avgPL_2 <- colMeans(PL_mat_2)
    
    # hold original distance of mean vectors and counter for test
    dist <- distance(avgPL_1, avgPL_2)
    count <- 0
    
    # merge the matrices and count number of rows in each
    merged <- rbind(PL_mat_1, PL_mat_2)
    n_1 <- nrow(PL_mat_1)
    n_2 <- nrow(PL_mat_2)
    total_rows <- n_1 + n_2
    
    # run k permutations
    for (i in 1:reps) {
      # create permutation of row indices
      permutation <- sample(1:total_rows)
      
      # get average PL vector for each new group
      avgPL_1 <- colMeans(merged[permutation[1:n_1],])
      avgPL_2 <- colMeans(merged[permutation[(n_1+1):total_rows],])
      
      # increase count if distance of this permutation is greater than original
      perm_dist <- distance(avgPL_1, avgPL_2)
      if (perm_dist > dist) {
        count <- count + 1
      }
    }
    
    # output p-value of comparison
    print(paste0("p-value for ", dir_1, " vs. ", dir_2, " -> ", count / reps))
  }
}


plot_PL <- function(PL, image_file, x_lim, y_lim, n, resolution) {
  # define color list for landscape functions
  colors <- c("#648FFF", "#6E77F8", "#785EF0", "#AA42B8", "#DC267F", "#ED4440", "#FE6100", "#FF8900", "#FFB000")
  
  # create plot
  p <- ggplot() +
    scale_x_continuous(expand=c(0, 0), limits=1.1*x_lim) + 
    scale_y_continuous(expand=c(0, 0), limits=1.1*y_lim) + 
    theme_classic() + xlab("(birth+death)/2") + ylab("(death-birth)/2")
  
  # plot each of the landscape functions
  for (i in 2:(n+1)) {
    temp <- data.frame(X=PL[,1], Y=PL[,i])
    p <- p + geom_line(data=temp, aes(X,Y), color=colors[i %% 8 + 1], linewidth=0.3)
  }
  
  # convert width to inches and write image file
  width <- resolution / 300    # 300 pixels per inch default
  ggsave(image_file, width=width, height=width*(y_lim[2]/x_lim[2]))
}


PLs_to_mat <- function(path, depth, all) {
  # create holder for vector representations of PLs
  PLs <- list()
  
  # get files
  files <- get_files(path, "\\PL.csv$")
  n <- length(files)
  
  # iterate through all PL files
  for (i in 1:n) {
    # read the CSV and add landscape to total
    PL_file <- concat_paths(path, files[[i]])
    PL <- data.matrix(read.csv(PL_file))[,-1]    # remove first column
    
    # determine how many functions to use
    m <- ncol(PL)
    if (all | m < depth) {
      PLs[[i]] <- as.vector(PL)    # use all discretized landscape functions
    } else {
        PLs[[i]] <- as.vector(PL[,1:depth])    # only use top 1 to depth functions
    }
  }
  
  # return a matrix with each row as a vector for each PL
  return(do.call(rbind, PLs))
}


get_files <- function(path, pattern) {
  # get list of discretized CSVs recursively
  files <- list.files(path, pattern=pattern, include.dirs=TRUE, recursive=TRUE)
}


make_dir <- function(path) {
  # check if directory exists and if not, make it
  if (!dir.exists(path)) {
    dir.create(path, recursive=TRUE)
  }
}


concat_paths <- function(path_1, path_2) {
  # make sure first path has file separator at the end
  n <- nchar(path_1)
  if (substr(path_1, n, n) != .Platform$file.sep) {
    path_1 <- paste0(path_1, .Platform$file.sep)
  }
  
  # return correct joined paths
  return(paste0(path_1, path_2))
}


distance <- function(u, v) {
  # return euclidean distance
  sqrt(sum((u - v) ^ 2))
}

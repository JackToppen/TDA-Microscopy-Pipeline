#Iryna Hartsock
#Pick the best patches based on cell number count
#March 2022

library(tictoc)

tic()
source("config.R")
#save_computations <- TRUE

# moved these into the config file
#save_file_location <- save_computation_file_location
#save_filename <- save_computation_filename

#compute number of cells for each file
cell_number_r <- list()
red_patch <- list()
red_patches <- list()
cell_number_gr <- list()
green_patch <- list()
green_patches <- list()
cell_number_pos <- list()
pos_patch <- list()
pos_patches <- list()
cell_number_neg <- list()
neg_patch <- list()
neg_patches <- list()
for (j in 1:ncohorts){                    #j=1 Gata6, j=2 HA
  cell_number_r[[j]] <- list()
  red_patch[[j]] <- list()
  red_patches[[j]] <- list()
  cell_number_gr[[j]] <- list()
  green_patch[[j]] <- list()
  green_patches[[j]] <- list()
  cell_number_pos[[j]] <- list()
  pos_patch[[j]] <- list()
  pos_patches[[j]] <- list()
  cell_number_neg[[j]] <- list()
  neg_patch[[j]] <- list()
  neg_patches[[j]] <- list()
  for (c in 1:nconcentrations){                  #c corresponds to 4 concentrations
    red_patch[[j]][[c]] <- list()
    cell_number_r[[j]][[c]]<- list()
    red_patches[[j]][[c]] <- list()
    green_patch[[j]][[c]] <- list()
    cell_number_gr[[j]][[c]]<- list()
    green_patches[[j]][[c]] <- list()
    pos_patch[[j]][[c]] <- list()
    cell_number_pos[[j]][[c]]<- list()
    pos_patches[[j]][[c]] <- list()
    neg_patch[[j]][[c]] <- list()
    cell_number_neg[[j]][[c]]<- list()
    neg_patches[[j]][[c]] <- list()
    for (i in 1:nsamples){               #k corresponds to 4 patches
      red_patch[[j]][[c]][[i]] <- list()
      cell_number_r[[j]][[c]][[i]] <- list()
      green_patch[[j]][[c]][[i]] <- list()
      cell_number_gr[[j]][[c]][[i]] <- list()
      pos_patch[[j]][[c]][[i]] <- list()
      cell_number_pos[[j]][[c]][[i]] <- list()
      neg_patch[[j]][[c]][[i]] <- list()
      cell_number_neg[[j]][[c]][[i]] <- list()
      for (k in 1:1){
        all_cells <- read.csv(file=paste0(data_file_location,data_filenames[[i+(c-1)*15+(j-1)*60]]), header=TRUE)
        all_cells <- all_cells[which( 500 <= all_cells[,1] & all_cells[,1] < 2500 & 500 <= all_cells[,2] & all_cells[,2] < 2500),]
       # if (k==1){
       #   all_cells <- all_cells[which( 500 <= all_cells[,1] & all_cells[,1] < 1500 & 500 <= all_cells[,2] & all_cells[,2] < 1500),]
       # } else if (k==2){
      #    all_cells <- all_cells[which( 1500 <= all_cells[,1] & all_cells[,1] < 2500 & 500 <= all_cells[,2] & all_cells[,2] < 1500),]
       # } else if (k==3){
        #  all_cells <- all_cells[which( 500 <= all_cells[,1] & all_cells[,1] < 1500 & 1500 <= all_cells[,2] & all_cells[,2] < 2500),]
       # } else {
       #   all_cells <- all_cells[which( 1500 <= all_cells[,1] & all_cells[,1] < 2500 & 1500 <= all_cells[,2] & all_cells[,2] < 2500),]
       # }
        #red
        red_patch[[j]][[c]][[i]][[k]] <- all_cells[ which((all_cells[,8]==0 & all_cells[,9]==1)), ] 
        cell_number_r[[j]][[c]][[i]][[k]] <- nrow(red_patch[[j]][[c]][[i]][[k]])
        #red_patches[[j]][[c]][[i+(k-1)*15]] <- red_patch[[j]][[c]][[k]][[i]] #list of all patches for a fixed concentration and cohort
        #green
        green_patch[[j]][[c]][[i]][[k]] <- all_cells[ which((all_cells[,8]==1 & all_cells[,9]==0)), ] 
        cell_number_gr[[j]][[c]][[i]][[k]] <- nrow(green_patch[[j]][[c]][[i]][[k]])
        #green_patches[[j]][[c]][[i+(k-1)*15]] <- green_patch[[j]][[c]][[k]][[i]]
        #double positive
        pos_patch[[j]][[c]][[i]][[k]] <- all_cells[ which((all_cells[,8]==1 & all_cells[,9]==1)), ] 
        cell_number_pos[[j]][[c]][[i]][[k]] <- nrow(pos_patch[[j]][[c]][[i]][[k]])
        #pos_patches[[j]][[c]][[i+(k-1)*15]] <- pos_patch[[j]][[c]][[k]][[i]] 
        #double negative
        neg_patch[[j]][[c]][[i]][[k]] <- all_cells[ which((all_cells[,8]==0 & all_cells[,9]==0)), ] 
        cell_number_neg[[j]][[c]][[i]][[k]] <- nrow(neg_patch[[j]][[c]][[i]][[k]])
        #neg_patches[[j]][[c]][[i+(k-1)*15]] <- neg_patch[[j]][[c]][[k]][[i]]  
      }
    }
  }
}


counts <- list()
for (j in 1:ncohorts){ 
  counts[[j]] <- list()
  for (c in 1:nconcentrations){
    counts[[j]][[c]] <- vector()
    for (i in 1:nsamples){
      counts[[j]][[c]] <- rbind(counts[[j]][[c]],c(cell_number_r[[j]][[c]][[i]][[1]], cell_number_gr[[j]][[c]][[i]][[1]], cell_number_pos[[j]][[c]][[i]][[k]], cell_number_neg[[j]][[c]][[i]][[k]]))
    }
  }
}

#Compute histograms
cell_number <- cell_number_r #pick a cell type
for (j in 1:2){
  if (j==1){
    cohort <- "Gata6"
  } else {cohort <- "HA"}
  for (c in 1:4){
    if (c==1){
      concentration <- 0
    } else if (c==2){
      concentration <- 5
    } else if (c==3){
      concentration <- 15
    } else {concentration <- 25 }
    
    h <- hist(unlist(cell_number[[j]][[c]]), breaks=20, ylim=c(0,35), main=sprintf("Histogram for %s concentration %1.0f, p.%1.0f", cohort, concentration, percentile), xlab="number of cells")
    text(h$mids, h$counts, labels=h$counts, adj=c(0.5,-0.5))
  }
}


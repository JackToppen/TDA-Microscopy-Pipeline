#Iryna Hartsock
#Pick the best patches based on cell number count

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
    for (i in 1:nimages){               #k corresponds to 4 patches
      all_cells <- read.csv(file=paste0(data_file_location,data_filenames[[i+(c-1)*nimages+(j-1)*nimages*nconcentrations]]), header=TRUE)
      red_patch[[j]][[c]][[i]] <- list()
      cell_number_r[[j]][[c]][[i]] <- list()
      green_patch[[j]][[c]][[i]] <- list()
      cell_number_gr[[j]][[c]][[i]] <- list()
      pos_patch[[j]][[c]][[i]] <- list()
      cell_number_pos[[j]][[c]][[i]] <- list()
      neg_patch[[j]][[c]][[i]] <- list()
      cell_number_neg[[j]][[c]][[i]] <- list()
      for (k in 1:sqrt(npatches)){
        for (r in 1:sqrt(npatches)){
          #every image is 3000x3000 
          patch <- all_cells[which( ((3000/sqrt(npatches))*(r-1)) <= all_cells[,1] 
                                    & all_cells[,1] < ((3000/sqrt(npatches))*r) 
                                    & ((3000/sqrt(npatches))*(k-1)) <= all_cells[,2] 
                                    & all_cells[,2] < ((3000/sqrt(npatches))*k)),]
          #red
          red_patch[[j]][[c]][[i]][[r+(k-1)*sqrt(npatches)]] <- patch[ which((patch[,8]==0 & patch[,9]==1)), ] 
          cell_number_r[[j]][[c]][[i]][[r+(k-1)*sqrt(npatches)]] <- nrow(red_patch[[j]][[c]][[i]][[r+(k-1)*sqrt(npatches)]])
          #red_patches[[j]][[c]][[i+(k-1)*15]] <- red_patch[[j]][[c]][[k]][[i]] #list of all patches for a fixed concentration and cohort
          #green
          green_patch[[j]][[c]][[i]][[r+(k-1)*sqrt(npatches)]] <- patch[ which((patch[,8]==1 & patch[,9]==0)), ] 
          cell_number_gr[[j]][[c]][[i]][[r+(k-1)*sqrt(npatches)]] <- nrow(green_patch[[j]][[c]][[i]][[r+(k-1)*sqrt(npatches)]])
          #green_patches[[j]][[c]][[i+(k-1)*15]] <- green_patch[[j]][[c]][[k]][[i]]
          #double positive
          pos_patch[[j]][[c]][[i]][[r+(k-1)*sqrt(npatches)]] <- patch[ which((patch[,8]==1 & patch[,9]==1)), ] 
          cell_number_pos[[j]][[c]][[i]][[r+(k-1)*sqrt(npatches)]] <- nrow(pos_patch[[j]][[c]][[i]][[r+(k-1)*sqrt(npatches)]])
          #pos_patches[[j]][[c]][[i+(k-1)*15]] <- pos_patch[[j]][[c]][[k]][[i]] 
          #double negative
          neg_patch[[j]][[c]][[i]][[r+(k-1)*sqrt(npatches)]] <- patch[ which((patch[,8]==0 & patch[,9]==0)), ] 
          cell_number_neg[[j]][[c]][[i]][[r+(k-1)*sqrt(npatches)]] <- nrow(neg_patch[[j]][[c]][[i]][[r+(k-1)*sqrt(npatches)]])
          #neg_patches[[j]][[c]][[i+(k-1)*15]] <- neg_patch[[j]][[c]][[k]][[i]]  
        }
      }
    }
  }
}
toc()


counts <- list()
for (j in 1:ncohorts){
  counts[[j]] <- cbind(matrix(unlist(cell_number_r[[j]]), nrow=length(unlist(cell_number_r[[j]]))), 
                       matrix(unlist(cell_number_gr[[j]]), nrow=length(unlist(cell_number_gr[[j]]))), 
                       matrix(unlist(cell_number_pos[[j]]), nrow=length(unlist(cell_number_pos[[j]]))), 
                       matrix(unlist(cell_number_neg[[j]]), nrow=length(unlist(cell_number_neg[[j]]))))
}

#Compute histograms
cell_number <- cell_number_r #pick a cell type
for (j in 1:ncohorts){
  if (j==1){
    cohort <- "Gata6"
  } else {cohort <- "HA"}
  for (c in 1:nconcentrations){
    if (c==1){
      concentration <- 0
    } else if (c==2){
      concentration <- 5
    } else if (c==3){
      concentration <- 15
    } else {concentration <- 25 }
    hist(unlist(cell_number[[j]][[c]]), breaks=40, ylim=c(0,300), main=sprintf("Histogram for %s concentration %1.0f, p.%1.0f", cohort, concentration, percentile), xlab="number of cells")
     #text(h$mids, h$counts, labels=h$counts, adj=c(0.5,-0.5))
  }
}

#Compute number of patches that are either empty or contain only one point
cell_number <- cell_number_r
almost_empty <- list()
for (j in 1:ncohorts){
  almost_empty[[j]] <- list()
  for (c in 1:nconcentrations){
    print(sprintf("concentration: %i", c))
    almost_empty[[j]][[c]] <- list()
    for (i in 1:nimages){
      almost_empty[[j]][[c]][[i]] <- vector()
      for (k in 1:npatches){
        if ((cell_number_r[[j]][[c]][[i]][[k]]+cell_number_gr[[j]][[c]][[i]][[k]]+cell_number_pos[[j]][[c]][[i]][[k]]+cell_number_neg[[j]][[c]][[i]][[k]]) <= 10){
          almost_empty[[j]][[c]][[i]] <- append(almost_empty[[j]][[c]][[i]], k)
        }
      }
      print(length(unlist(almost_empty[[j]][[c]][[i]])))
    }
  }
}

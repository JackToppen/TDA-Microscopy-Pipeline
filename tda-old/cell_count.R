#Iryna Hartsock
#Pick the best patches based on cell number count

library(tictoc)

tic()
source("config.R")
#save_computations <- TRUE

# moved these into the config file
#save_file_location <- save_computation_file_location
#save_filename <- save_computation_filename

#compute number of cells for every patch
cell_number_r <- list()
cell_number_gr <- list()
cell_number_pos <- list()
cell_number_neg <- list()
for (j in 1:ncohorts){                    #j=1 Gata6, j=2 HA
  cell_number_r[[j]] <- list()
  cell_number_gr[[j]] <- list()
  cell_number_pos[[j]] <- list()
  cell_number_neg[[j]] <- list()
  for (c in 1:nconcentrations){                  #c corresponds to 4 concentrations
    cell_number_r[[j]][[c]]<- list()
    cell_number_gr[[j]][[c]]<- list()
    cell_number_pos[[j]][[c]]<- list()
    cell_number_neg[[j]][[c]]<- list()
    for (i in 1:nimages){               #k corresponds to 4 patches
      all_cells <- read.csv(file=paste0(data_file_location,data_filenames[[i+(c-1)*nimages+(j-1)*nimages*nconcentrations]]), header=TRUE)
      cell_number_r[[j]][[c]][[i]] <- list()
      cell_number_gr[[j]][[c]][[i]] <- list()
      cell_number_pos[[j]][[c]][[i]] <- list()
      cell_number_neg[[j]][[c]][[i]] <- list()
      for (k in 1:sqrt(npatches)){
        for (r in 1:sqrt(npatches)){
          #every image is 2850x2850 
          patch <- all_cells[which( ((2850/sqrt(npatches))*(r-1)) <= all_cells[,1] 
                                    & all_cells[,1] < ((2850/sqrt(npatches))*r) 
                                    & ((2850/sqrt(npatches))*(k-1)) <= all_cells[,2] 
                                    & all_cells[,2] < ((2850/sqrt(npatches))*k)),]
          #red
          red_patch <- patch[ which((patch[,8]==1)), ] 
          cell_number_r[[j]][[c]][[i]][[r+(k-1)*sqrt(npatches)]] <- nrow(red_patch)
          #green
          green_patch <- patch[ which((patch[,8]==10)), ] 
          cell_number_gr[[j]][[c]][[i]][[r+(k-1)*sqrt(npatches)]] <- nrow(green_patch)
          #double positive
          pos_patch <- patch[ which((patch[,8]==11)), ] 
          cell_number_pos[[j]][[c]][[i]][[r+(k-1)*sqrt(npatches)]] <- nrow(pos_patch)
          #double negative
          neg_patch <- patch[ which((patch[,8]==0)), ] 
          cell_number_neg[[j]][[c]][[i]][[r+(k-1)*sqrt(npatches)]] <- nrow(neg_patch)
        }
      }
    }
  }
}
toc() 


q_rG <- quantile(unlist(cell_number_r[[1]]), prob=seq(0,1,1))
q_rH <- quantile(unlist(cell_number_r[[2]]), prob=seq(0,1,1))
q_r <- quantile(unlist(cell_number_r), prob=seq(0,1,0.1))
q_gr <- quantile(unlist(cell_number_gr), prob=seq(0,1,0.1))
q_pos <- quantile(unlist(cell_number_pos), prob=seq(0,1,0.1))
q_neg <- quantile(unlist(cell_number_neg), prob=seq(0,1,0.1))

l <- unlist(cell_number_gr[[1]][[4]])
sum(l)
print(length(l[which(l<=1000)]))


counts <- list()
for (j in 1:ncohorts){
  counts[[j]] <- cbind(matrix(unlist(cell_number_r[[j]]), nrow=length(unlist(cell_number_r[[j]]))), 
                       matrix(unlist(cell_number_gr[[j]]), nrow=length(unlist(cell_number_gr[[j]]))), 
                       matrix(unlist(cell_number_pos[[j]]), nrow=length(unlist(cell_number_pos[[j]]))), 
                       matrix(unlist(cell_number_neg[[j]]), nrow=length(unlist(cell_number_neg[[j]]))))
}


#compute number of cells for every image
cell_number_r <- list()
cell_number_gr <- list()
cell_number_pos <- list()
cell_number_neg <- list()
intensities_r <- list()
intensities_gr <- list()
red <- list()
for (j in 1:ncohorts){                    #j=1 Gata6, j=2 HA
  cell_number_r[[j]] <- list()
  cell_number_gr[[j]] <- list()
  cell_number_pos[[j]] <- list()
  cell_number_neg[[j]] <- list()
  intensities_r[[j]] <- list()
  intensities_gr[[j]] <- list()
  red[[j]] <- list()
  for (i in 1:(length(data_filenames)/2)){               #k corresponds to 4 patches
    all_cells <- read.csv(file=paste0(data_file_location,data_filenames[[i+(j-1)*60]]), header=TRUE)
    print(i+(j-1)*60)
    #red
    intensities_r[[j]][[i]] <- all_cells[,7]
    red[[j]][[i]] <- nrow(all_cells[ which(all_cells[,7] >= 0.42725448 & all_cells[,6] < 0.22028941), ] )
    red_patch <- all_cells[ which((all_cells[,8]==0 & all_cells[,9]==1)), ] 
    cell_number_r[[j]][[i]] <- nrow(red_patch)
    #green
    intensities_gr[[j]][[i]] <- all_cells[,6]
    green_patch <- all_cells[ which((all_cells[,8]==1 & all_cells[,9]==0)), ] 
    cell_number_gr[[j]][[i]] <- nrow(green_patch)
    #double positive
    pos_patch <- all_cells[ which((all_cells[,8]==1 & all_cells[,9]==1)), ] 
    cell_number_pos[[j]][[i]] <- nrow(pos_patch)
    #double negative
    neg_patch <- all_cells[ which((all_cells[,8]==0 & all_cells[,9]==0)), ] 
    cell_number_neg[[j]][[i]] <- nrow(neg_patch)
  }
}

sum(unlist(red[[1]][46:60]))
q_r <- quantile(unlist(intensities_r[[1]][1:15] ), prob=seq(0,1,0.05))
q_gr <- quantile(unlist(intensities_gr[[1]][46:60] ), prob=seq(0,1,0.05))

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
    hist((unlist(cell_number_r[[j]][[c]])+unlist(cell_number_pos[[j]][[c]])), breaks=20, ylim=c(0,50), main=sprintf("Histogram for %s concentration %1.0f, p.%1.0f", cohort, concentration, percentile), xlab="number of cells")
     #text(h$mids, h$counts, labels=h$counts, adj=c(0.5,-0.5))
  }
}



par(pty="s")
plot(density(unlist(cell_number_r[[1]][[1]]), bw=100), col= "red", xlab="number of cells", xlim=c(0,4000))
lines(density(unlist(cell_number_r[[1]][[2]]), bw=100), col="purple")
lines(density(unlist(cell_number_r[[1]][[3]]), bw=100), col="blue")
lines(density(unlist(cell_number_r[[1]][[4]]), bw=100), col="green")
legend("topright", c("0","5", "15", "25"), lty = c(1,1), col = c("red","purple", "blue", "green"))


counts <- list(cell_number_r, cell_number_gr, cell_number_pos, cell_number_neg)
y <- list(0.0021, 0.0025, 0.0007, 0.0006)
for (i in 1:4){
  cell_number <- counts[[i]]
  max_element <- list()
  data <- list()
  for (j in 1:2){
    data[[j]] <- data.frame(
      concentration = c(rep("0", length(unlist(cell_number[[j]][[1]]))), rep("5", length(unlist(cell_number[[j]][[2]]))), rep("15", length(unlist(cell_number[[j]][[3]]))), rep("25", length(unlist(cell_number[[j]][[4]])))),
      cell_count = c( unlist(cell_number[[j]][[1]]), unlist(cell_number[[j]][[2]]), unlist(cell_number[[j]][[3]]), unlist(cell_number[[j]][[4]]))
    )
    max_element[[j]] <- max(c( unlist(cell_number[[j]][[1]]), unlist(cell_number[[j]][[2]]), unlist(cell_number[[j]][[3]]), unlist(cell_number[[j]][[4]])))
  }
  max_count <- max(max_element[[1]], max_element[[2]])
  for (j in 1:2){
    ggplot(data[[j]], aes(x = cell_count, fill = concentration)) + geom_density(alpha = 0.4)+scale_fill_discrete(breaks=c('0', '5', '15', '25'))+coord_cartesian(ylim = c(0, y[[i]]), xlim=c(0,max_count))
  }
}




# Represent it
  ggplot( aes(x=value, fill=type)) +
  geom_density( color="#e9ecef", alpha=0.6) +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_ipsum() +
  labs(fill="")

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

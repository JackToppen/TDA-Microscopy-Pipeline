#Ashleigh Thomas, Iryna Hartsock
# function to plot landscapes differently from tda-tools
# input can contain ...
# new way of plotting with specified maximum y value and no labels
# code edited from original tda-tools code by Jose Bouza

plot_landscape_internal <- function(internal, y_max, y_min=0){
  infinity_sub=-1
  line_width <- 1
  # colors <- colorRamp(c("red", "blue"), numLevels(internal)) # outputs a function
  # colors <- myColorRamp(c("red", "blue"), numLevels(internal)) # outputs a vector
  colors <- c(1:numLevels(internal)) # outputs a vector
  
  for(level in 1:numLevels(internal)){
    level_d <- accessLevel(internal,level)
    if(infinity_sub != -1){
      level_d[level_d == Inf] <- infinity_sub
      level_d[level_d == -Inf] <- -infinity_sub
    }
  }
  level1 <- accessLevel(internal,1)
  
  # plot first level
  if (missing(y_max)){ # plot as usual
    plot(level1[,1],level1[,2], type='l', xlab='', ylab='', bty="n", col=plasma(1), lwd=line_width)
  } else { # specify limits of y axis
    plot(level1[,1],level1[,2], type='l', ann=FALSE, bty="n", col=plasma(1), lwd=line_width, ylim=c(y_min,y_max))#viridis(numLevels(internal))
  }
  
  if(numLevels(internal) > 1){
    for(level in 2:numLevels(internal)){
      level_d <- accessLevel(internal,level)
      lines(level_d[,1], level_d[,2], col=plasma(15)[level %% 15], lwd=line_width)
    }
  }
}

plot_landscape <- function(PersistenceLandscape, y_max){
  internal <- PersistenceLandscape$getInternal()
  plot_landscape_internal(internal, y_max)
}


plot_landscape_internal_levels <- function(internal, l1, l2, y_max, y_min=0){
  infinity_sub=-1
  line_width <- 1
   mycolors <- rep(c("midnightblue","slateblue", "royalblue2", "deepskyblue4", "deepskyblue3",
                         "turquoise4", "chartreuse4","olivedrab4", "olivedrab3", "yellowgreen", 
                     "yellow3","gold3", "gold2","gold1", "gold" ), 2) # outputs a function
  # colors <- myColorRamp(c("red", "blue"), numLevels(internal)) # outputs a vector
  #colors <- c(l1:l2) # outputs a vector
    #"gold3", "goldenrod", "tan2", "tan1"
  for(level in l1:l2){
    level_d <- accessLevel(internal,level)
    if(infinity_sub != -1){
      level_d[level_d == Inf] <- infinity_sub
      level_d[level_d == -Inf] <- -infinity_sub
    }
  }
  level1 <- accessLevel(internal,l1)
  
  # plot first level
  if (missing(y_max)){ # plot as usual
    plot(level1[,1],level1[,2], type='l', xlab='', ylab='',  col=mycolors[1], lwd=line_width)
  } else { # specify limits of y axis
    plot(level1[,1],level1[,2], type='l', ann=FALSE, col=mycolors[1], lwd=line_width, ylim=c(y_min,y_max))#viridis(numLevels(internal))
  }
  
  if( l2>l1){
    for(level in (l1+1):l2){
      level_d <- accessLevel(internal,level)
      lines(level_d[,1], level_d[,2], col=mycolors[level], lwd=line_width)
    }
    
  }
}

#plot certain landscapes' levels
plot_landscape_levels <- function(PersistenceLandscape, l1, l2, y_max){
  internal <- PersistenceLandscape$getInternal()
  plot_landscape_internal_levels(internal, l1, l2, y_max)
}



get_landscape_y_max <- function(PersistenceLandscape){
  y_max <- 0
  internal <- PersistenceLandscape$getInternal()
  level1 <- accessLevel(internal,1)
  if (length(level1)==2) {
    y_max <- max(y_max,level1[2])
  } else {
    y_max <- max(level1[,2])
  }
  return(y_max)
}


vectorize_landscapes <- function(PL_list, depth_cap=0){
  # input: list of landscapes, highest depth of landscape to include in vector
  # landscape data is a collection of (x,y) points that approximate 
  # the persistence landscape of some data. 
  # vectorize_landscapes returns a matrix where each row is the concatenation 
  # of the y values for each level (up to depth_cap) of a landscape. 
  if (depth_cap == 0){ # no depth cap -- set to max depth
    max_depth <- 0
    for (i in 1:length(PL_list)){
      max_depth <- max(max_depth, dim(PL_list[[i]]$getInternal())[1])
    }
    depth_cap <- max_depth
  }
  
  vect_length <- depth_cap*dim(PL_list[[1]]$getInternal())[2]
  vect_PLs <- matrix(0, nrow=length(PL_list), ncol=vect_length)
  for (i in 1:length(PL_list)){
    if (dim(PL_list[[i]]$getInternal())[1] < depth_cap) { # fewer than depth_cap levels in this landscape
      # need the transpose because as.vector takes columns of a matrix, not rows
      temp_vec <- as.vector(t(PL_list[[i]]$getInternal()[,,2]))
      vect_PLs[i,1:length(temp_vec)] <- temp_vec
    } else {
      temp_vec <- as.vector(t(PL_list[[i]]$getInternal()[1:depth_cap,,2]))
      vect_PLs[i,] <- temp_vec
    }
  }
  return(vect_PLs)
}


# computes landscape from persistence diagram
# returns 0 landscape if persistence diagram is empty (instead of erroring)
landscape0 <- function(data, degree, exact=FALSE, dx, min_x, max_x){
  if (length(data)==0) { # empty persistence diagram
    tdatools::landscape(matrix(0, nrow=1,ncol=2), degree=degree, exact=exact, dx=dx, min_x=min_x, max_x=max_x)
  } else {
    tdatools::landscape(data, degree=degree, exact=exact, dx=dx, min_x=min_x, max_x=max_x)
  }
}

# functions to convert persistence diagrams between data structures
# GUDHI (Dionysus) --> tda-tools (Ripser)
gudhi2tdatools <- function(gudhi_dgm) {
  # extract first homology from gudhi diagram
  
  if (length(which(gudhi_dgm[,1]==1)) >1){
    h1 <- gudhi_dgm[which(gudhi_dgm[,1]==1),][,2:3] # size is dim(h1) x 2
     #plot(h1, asp=1)
  }
  if (length(which(gudhi_dgm[,1]==1))==1){
    h <- gudhi_dgm[which(gudhi_dgm[,1]==1),2:3] # size is dim(h1) x 2
    h1 <- rbind(rep(0,2),h)
     #plot(h1, asp=1)
  }
  if (length(which(gudhi_dgm[,1]==1))==0){
    h1 <- rbind(rep(0,2),rep(0,2))
  }
  # extract zero homology from gudhi diagram
  if (length(which(gudhi_dgm[,1]==0)) >1){
    h0 <- gudhi_dgm[which(gudhi_dgm[,1]==0),][,2:3] # size is dim(h1) x 2
     #plot(h0, asp=1)
  }
  if (length(which(gudhi_dgm[,1]==0))==1){
    h <- gudhi_dgm[which(gudhi_dgm[,1]==0),2:3] # size is dim(h1) x 2
    h0 <- rbind(rep(0,2),h)
    # plot(h1, asp=1)
  }
  if (length(which(gudhi_dgm[,1]==0))==0){
    h0 <- rbind(rep(0,2),rep(0,2))
  }
  # create tdatools diagram; load with zero and first homology
  tdatools_dgm <- diagram(0,'point-cloud', dim_max = 1)
  tdatools_dgm$pairs[[2]] <- h1
  tdatools_dgm$pairs[[1]] <- h0
  
  return(tdatools_dgm)
}

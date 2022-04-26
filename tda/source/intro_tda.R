### Part from the Introduction to Topological Data Analysis with R 
### Peter Bubenik
### October 11, 2019

#install.packages("TDA") # R Topological Data Analysis package
#install.packages("deldir") # R Delaunay and Voronoi package
#install.packages("kernlab") # R Support Vector Machine package
library(TDA) 
library(deldir) 
library("Matrix") # package for sparse matrices
library("kernlab")
par(pty="s") # force the plotting region to be square

###############################
########## Functions ##########
###############################

# Euclidean distance between vectors
euclidean.distance <- function(u, v) sqrt(sum((u - v) ^ 2))

# Plot the Voronoi cells and dual and Delaunay complex
plot.delaunay <- function(X){
  DelVor <- deldir(X[,1], X[,2], suppressMsge = TRUE)
  # Voronoi cells:
  plot(DelVor, pch=20, col=c("black","red","blue"), wlines= ("tess"), asp=1)
  # Voronoi cells and their dual (the Delaunay complex):
  plot(DelVor, pch=20, col=c("black","red","blue"), wlines= ("both"), asp=1)
  # Delaunay complex:
  plot(DelVor, pch=20, col=c("black","red","blue"), wlines= ("triang"), asp=1)
}

# Plot representative cycles for Delaunay complex
plot.delaunay.cycle <- function(X){
  PH.output <- alphaComplexDiag(X, maxdimension = 1, library = c("GUDHI", "Dionysus"), 
                                location = TRUE)
  PD <- PH.output[["diagram"]]
  ones <- which(PD[, 1] == 1)
  persistence <- PD[ones,3] - PD[ones,2]
  cycles <- PH.output[["cycleLocation"]][ones[order(persistence)]]
  for (i in 1:length(cycles)){
    plot(X, pch=20, col="blue", asp=1)
    for (j in 1:dim(cycles[[i]])[1])
      lines(cycles[[i]][j,,])
  }
}


# Plot Persistence Landscape
plot.landscape <- function(PL,t.vals){
  plot(t.vals,PL[1,],type="l",ylab="Persistence",xlab="Parameter values",col=1,ylim=c(min(PL),max(PL)))
  for(i in 2:dim(PL)[1])
    lines(t.vals,PL[i,],type="l",col=i)
}

# Matrix of persistence landscape row vectors from list of persistence landscapes
landscape.matrix.from.list <- function(PL.list){
  n <- length(PL.list)
  m <- ncol(PL.list[[1]])
  max.depth <- integer(n)
  for (i in 1:n)
    max.depth[i] <- nrow(PL.list[[i]])
  K <- max(max.depth)
  PL.matrix <- Matrix(0, nrow = n, ncol = m*K, sparse = TRUE)
  for (i in 1:n)
    for (j in 1:max.depth[i])
      PL.matrix[i,(1+(j-1)*m):(j*m)] <- PL.list[[i]][j,]
  return(PL.matrix)
}

# Convert a vector to a persistence landscape
landscape.from.vector <- function(PL.vector, t.vals){
  m <- length(t.vals)
  K <- length(PL.vector)/m
  PL <- Matrix(0, nrow = K, ncol=m, sparse = TRUE)
  for (i in 1:K){
    PL[i,1:m] <- PL.vector[(1+(i-1)*m):(i*m)]
  }
  return(PL)
}

# Take difference of vectors of potentially different lengths
difference.vectors <-function(vector.1,vector.2){
  length.1 <- length(vector.1)
  length.2 <- length(vector.2)
  difference.vector = numeric(max(length.1,length.2))
  difference.vector = as(difference.vector, "sparseVector")
  difference.vector[1:length.1] = difference.vector[1:length.1] + vector.1
  difference.vector[1:length.2] = difference.vector[1:length.2] - vector.2
}

# Permutation test for two matrices consisting of row vectors
permutation.test <- function(M1 ,M2, num.repeats = 10000){
  # append zeros if necessary so that the matrices have the same number of columns
  num.columns <- max(ncol(M1),ncol(M2))
  M1 <- cbind(M1, Matrix(0,nrow=nrow(M1),ncol=num.columns-ncol(M1)))
  M2 <- cbind(M2, Matrix(0,nrow=nrow(M2),ncol=num.columns-ncol(M2)))
  t.obs <- euclidean.distance(colMeans(M1),colMeans(M2))
  k <- dim(M1)[1]
  M <- rbind(M1,M2)
  n <- dim(M)[1]
  count <- 0
  for (i in 1:num.repeats){
    permutation <- sample(1:n)
    t <- euclidean.distance(colMeans(M[permutation[1:k],]),colMeans(M[permutation[(k+1):n],]))
    if (t >= t.obs)
      count <- count + 1
  }
  return(count/num.repeats)
}

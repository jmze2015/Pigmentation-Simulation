####################################################
## ----2D Ising Model Pigmentation Application----##
## -------------(Jose Muñoz-Lopez)----------------##
####################################################

## Lattice sizes chosen N = 10,20,30,40 (or anything you want)
## Main Loop for N in the given vector

## Parameters

h <- 0
J1 <- 1
J2 <- -0.1
R1 <- 1
R2 <- 6
N <- 60


## Define a lattice with random spins (+1, -1) of size N^2
lattice <- matrix(sample(c(-1,1), N^2, replace = TRUE),ncol = N)

## Update Algorithm
Update <- function(A, h = 0, J1 = 1, J2 = -0.1, R1 = 1, R2 = 6){
  ## Pick site
  i <- sample(seq(1, nrow(A)),1)
  j <- sample(seq(1, ncol(A)),1)
  
  ## Squared distances
  d2 <- (row(A)-i)^2 + (col(A)-j)^2
  
  ## J1 term computation
  
  mask_j1 <- (d2 > 0) & (d2 <= R1^2)
  Term1 <- J1 * sum(A[mask_j1])
  
  ## J2 term computation
  
  mask_j2 <- (d2 > R1^2) & (d2 < R2^2)
  Term2 <- J2 * sum(A[mask_j2])

  
  A[i,j] <- sign(h + Term1 + Term2)
  
  return(A)
}



Sweep <- function(A,l){
  k <- N*N*l
  for (i in 1:k){
    A <- Update(A)
  }
  return(A)
}


## Visualization of a config.
Viz <- function(A){
  cols<- c("white", "black")
  image(
    1:nrow(A), 1:nrow(A),
    t(apply(A, 2, rev)),
    col = cols,
    breaks = c(-1, 0, 1),
    xlab = "",
    ylab = ""
  )  
}














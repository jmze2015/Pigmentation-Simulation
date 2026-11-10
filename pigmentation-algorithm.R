####################################################
## ----2D Ising Model Pigmentation Application----##
## -------------(Jose Muñoz-Lopez)----------------##
####################################################

## Lattice sizes chosen N = 60
## Main Loop for N in the given vector

## Parameters from Bio paper

# h <- 0
# J1 <- 1
# J2 <- -0.1
# R1 <- 1
# R2 <- 6

## Square Lattice Dimension
N <- 60
## Rectangular Lattice Dimension
M <- 80


## Define a lattice with random spins (+1, -1) of size N^2
lattice <- matrix(sample(c(-1,1), N^2, replace = TRUE),ncol = N)

## Define a rectangular CA (M x N)
rect <- matrix(sample(c(1,-1), N*M, replace = TRUE), ncol = N, nrow = M)

## Update Algorithm, regular boundary not periodic
Update <- function(A, h = 6, J1 = 1, J2 = -0.1, R1 = 1, R2 = 6){
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


## Visualization of a config.
Viz <- function(A, h, R1, R2){
  cols <- c("white", "black")
  
  image(
    1:ncol(A), 1:nrow(A),
    t(apply(A, 2, rev)),
    col = cols,
    breaks = c(-1, 0, 1),
    xlab = "",
    ylab = "",
    main = paste("(h =", h, ", R1 =", R1,
                 ", R2 =", R2, ")"),
    asp = 1
  )
}


pigment_simulation <- function(A, h = 6, J1 = 1, J2 = -0.1, R1 = 1, R2 = 6, total_updates, save_every){
  
  dir.create("frames", showWarnings = FALSE)
  
  unlink(list.files(
    "frames",
    pattern = "\\.png$",
    full.names = TRUE
  ))
  
  for (i in 1:total_updates){
    A <- Update(A, h, J1, J2, R1, R2)
    
    
    if (i %% save_every ==0){
      png(
        filename = sprintf("frames/frame_%08d.png", i),
        width = 800,
        height = 800
      )
      
      Viz(A, h, R1, R2)
      
      dev.off()
    }
  } 
  frames <- sort(list.files(
    "frames",
    pattern = "\\.png$",
    full.names = TRUE
  ))
  
  av_encode_video(
    input = frames,
    output = sprintf("h-%g_R1-%g_R2-%g_simulation_every%g.mp4", h, R1, R2, save_every),
    framerate = 30
  )
  
  unlink(frames)
}

## pigment_simulation(lattice, h = -6, R1 = 1, total_updates = 8*N*N, save_every = 200)
pigment_simulation(lattice, h = -3, R1 = 1, total_updates = 8*N*N, save_every = 200)
pigment_simulation(lattice, h = -1, R1 = 1, total_updates = 8*N*N, save_every = 200)
pigment_simulation(lattice, h = 0, R1 = 1, total_updates = 8*N*N, save_every = 200)
pigment_simulation(lattice, h = 1, R1 = 1, total_updates = 8*N*N, save_every = 200)
pigment_simulation(lattice, h = 3, R1 = 1, total_updates = 8*N*N, save_every = 200)
pigment_simulation(lattice, h = 6, R1 = 1, total_updates = 8*N*N, save_every = 200)




### Same ideas below but instead using a periodic boundary ###

## Update Algorithm, periodic boundary conditions
Update_periodic <- function(A, h = 6, J1 = 1, J2 = -0.1,
                            R1 = 1, R2 = 6){
  ## Lattice dimensions
  M <- nrow(A)
  N <- ncol(A)
  
  ## Pick site
  i <- sample(seq_len(M), 1)
  j <- sample(seq_len(N), 1)
  
  ## Periodic distances
  dr <- abs(row(A) - i)
  dc <- abs(col(A) - j)
  
  dr <- pmin(dr, M - dr)
  dc <- pmin(dc, N - dc)
  
  ## Squared periodic distances
  d2 <- dr^2 + dc^2
  
  ## J1 term computation
  mask_j1 <- (d2 > 0) & (d2 <= R1^2)
  Term1 <- J1 * sum(A[mask_j1])
  
  ## J2 term computation
  mask_j2 <- (d2 > R1^2) & (d2 < R2^2)
  Term2 <- J2 * sum(A[mask_j2])
  
  ## Update spin
  A[i,j] <- sign(h + Term1 + Term2)
  
  return(A)
}

## Visualization of a config.
Viz <- function(A, h, R1, R2){
  cols <- c("white", "black")
  
  image(
    1:ncol(A), 1:nrow(A),
    t(apply(A, 2, rev)),
    col = cols,
    breaks = c(-1, 0, 1),
    xlab = "",
    ylab = "",
    main = paste("Periodic Boundary (h =", h, ", R1 =", R1,
                 ", R2 =", R2, ")"),
    asp = 1
  )
}


periodic_pigment_simulation <- function(A, h = 6, J1 = 1, J2 = -0.1, R1 = 1, R2 = 6, total_updates, save_every){
  
  dir.create("frames", showWarnings = FALSE)
  
  unlink(list.files(
    "frames",
    pattern = "\\.png$",
    full.names = TRUE
  ))
  
  for (i in 1:total_updates){
    A <- Update_periodic(A, h, J1, J2, R1, R2)
    
    if (i %% save_every ==0){
      png(
        filename = sprintf("frames/frame_%08d.png", i),
        width = 800,
        height = 800
      )
      
      Viz(A, h, R1, R2)
      
      dev.off()
    }
  } 
  frames <- sort(list.files(
    "frames",
    pattern = "\\.png$",
    full.names = TRUE
  ))
  
  av_encode_video(
    input = frames,
    output = sprintf("h-%g_R1-%g_R2-%g_periodic_simulation_every%g.mp4", h, R1, R2, save_every),
    framerate = 30
  )
  unlink(frames)
}

periodic_pigment_simulation(lattice, h = 0, R1 = 1, total_updates = 8*N*N, save_every = 200)

## no real difference noted ##

### Below we fit 6 animations simulataneously ###

## Visualization of a configuration
Viz <- function(A, h) {
  cols <- c("white", "black")
  
  image(
    1:ncol(A), 1:nrow(A),
    t(apply(A, 2, rev)),
    col = cols,
    breaks = c(-1, 0, 1),
    xlab = "",
    ylab = "",
    main = paste("h =", h),
    asp = 1
  )
}

## Six-panel pigment simulation
pigment_simulation_6 <- function(A, J1 = 1, J2 = -0.1,
                                 R1 = 1, R2 = 6,
                                 total_updates, save_every) {
  
  ## Values of h
  h_vals <- c(-6, -3, 0, 1, 3, 6)
  
  ## Initialize six copies of the same lattice
  lattices <- replicate(6, A, simplify = FALSE)
  
  ## Create and clear frames folder
  dir.create("frames", showWarnings = FALSE)
  
  unlink(list.files(
    "frames",
    pattern = "\\.png$",
    full.names = TRUE
  ))
  
  ## Main simulation loop
  for (i in 1:total_updates) {
    
    ## Update each lattice
    for (k in 1:6) {
      lattices[[k]] <- Update(
        lattices[[k]],
        h = h_vals[k],
        J1 = J1,
        J2 = J2,
        R1 = R1,
        R2 = R2
      )
    }
    
    ## Save combined frame
    if (i %% save_every == 0) {
      
      png(
        filename = sprintf("frames/frame_%08d.png", i),
        width = 1200,
        height = 1600
      )
      
      par(mfrow = c(3, 2), mar = c(1, 1, 3, 1))
      
      for (k in 1:6) {
        Viz(lattices[[k]], h_vals[k])
      }
      
      dev.off()
    }
  }
  
  ## Collect frames
  frames <- sort(list.files(
    "frames",
    pattern = "\\.png$",
    full.names = TRUE
  ))
  
  ## Create MP4
  av::av_encode_video(
    input = frames,
    output = sprintf(
      "six_panel_R1-%g_R2-%g_every%g.mp4",
      R1, R2, save_every
    ),
    framerate = 30
  )
  
  ## Delete PNG frames
  unlink(frames)
  
  return(lattices)
}

pigment_simulation_6(
  lattice,
  R1 = 1,
  R2 = 6,
  total_updates = 8*N*N,
  save_every = 200
)



########## Visuals for the presentation Only ###########


## Visualization of a config.
Pattern_Viz <- function(A){
  nr <- nrow(A)
  nc <- ncol(A)
  
  cols <- c("white", "grey20", "#E5484D", "#4285F4")
  
  image(
    1:nc, 1:nr,
    t(apply(A, 2, rev)),
    col = cols,
    breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5),
    xlim = c(0.5, nc + 0.5),
    ylim = c(0.5, nr + 0.5),
    xlab = "",
    ylab = "",
    main = "Visualizing Discrete Neighborhoods",
    asp = 1,
    axes = FALSE,
    xaxs = "i",
    yaxs = "i"
  )
  
  ## Vertical grid lines
  segments(
    x0 = seq(0.5, nc + 0.5, 1),
    y0 = 0.5,
    x1 = seq(0.5, nc + 0.5, 1),
    y1 = nr + 0.5,
    col = "gray80",
    lwd = 0.5
  )
  
  ## Horizontal grid lines
  segments(
    x0 = 0.5,
    y0 = seq(0.5, nr + 0.5, 1),
    x1 = nc + 0.5,
    y1 = seq(0.5, nr + 0.5, 1),
    col = "gray80",
    lwd = 0.5
  )
  
  ## Outer border aligned with the lattice
  rect(0.5, 0.5, nc + 0.5, nr + 0.5, border = "black")
}

## Define matrix for visualization
A <- matrix(0, 20,20)

R1 <- 1
R2 <- 6

## Pick site
i <- sample(seq(1, nrow(A)),1)
j <- sample(seq(1, ncol(A)),1)

A[i,j] <- 1

## Squared distances
d2 <- (row(A)-i)^2 + (col(A)-j)^2

## J1 term computation

mask_j1 <- (d2 > 0) & (d2 <= R1^2)
A[mask_j1] <- 2
## J2 term computation

mask_j2 <- (d2 > R1^2) & (d2 < R2^2)
A[mask_j2] <- 3

Pattern_Viz(A)



#### Periodic Boundary Condition ##########

## Visualization of periodic neighborhoods
Pattern_Viz_Periodic <- function(N = 20, M = 20, R1 = 1, R2 = 6,
                        i = sample(1:N, 1),
                        j = sample(1:M, 1)) {
  
  ## Define lattice
  A <- matrix(0, N, M)
  
  ## Periodic distances
  dx <- abs(row(A) - i)
  dy <- abs(col(A) - j)
  
  dx <- pmin(dx, N - dx)
  dy <- pmin(dy, M - dy)
  
  d2 <- dx^2 + dy^2
  
  ## J1 neighborhood (red)
  mask_j1 <- (d2 > 0) & (d2 <= R1^2)
  A[mask_j1] <- 2
  
  ## J2 neighborhood (blue)
  mask_j2 <- (d2 > R1^2) & (d2 < R2^2)
  A[mask_j2] <- 3
  
  ## Selected site (black)
  A[i, j] <- 1
  
  ## Colors
  cols <- c("white", "grey20", "#E5484D", "#4285F4")
  
  ## Plot lattice
  image(
    1:M, 1:N,
    t(apply(A, 2, rev)),
    col = cols,
    breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5),
    xlim = c(0.5, M + 0.5),
    ylim = c(0.5, N + 0.5),
    xlab = "",
    ylab = "",
    main = "Periodic Neighborhoods",
    asp = 1,
    axes = FALSE,
    xaxs = "i",
    yaxs = "i"
  )
  
  ## Vertical grid lines
  segments(
    x0 = seq(0.5, M + 0.5, 1),
    y0 = 0.5,
    x1 = seq(0.5, M + 0.5, 1),
    y1 = N + 0.5,
    col = "gray80",
    lwd = 0.5
  )
  
  ## Horizontal grid lines
  segments(
    x0 = 0.5,
    y0 = seq(0.5, N + 0.5, 1),
    x1 = M + 0.5,
    y1 = seq(0.5, N + 0.5, 1),
    col = "gray80",
    lwd = 0.5
  )
  
  ## Outer border
  rect(0.5, 0.5, M + 0.5, N + 0.5, border = "black")
  
  invisible(A)
}


Pattern_Viz_Periodic(i = 2, j = 2)









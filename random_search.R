rm(list=ls()); graphics.off(); cat("\014")

cat("> Importing R packages..."); library(transport); library(DiceKriging); library(lhs); library(mvtnorm); library(reticulate); cat("done!\n")
cat("> Importing Python functions..."); source_python("surrogate_simulator.py"); cat("done!\n")

cat("\n> *** Importing utility for permutation-invariant objective1..."); source("utility.R"); cat("done!\n")


#*********************************************************************************
# Experimental setting
#*********************************************************************************

# ***************************************************
# Problem setup
# ***************************************************

m <- 5        # number of turbines
rho <- 0.1512 # from 2 * rotor diam / farm length = (2*126)/(333.33*5)

nfereps <- 5            # repetitions of function evaluation for each query
nSeeds <- 5             # just like in BO experiments
N <- 300                # just like in BO experiments

ref_mean <- c(-1,-1)      # bi-variate mean of the reference point-cloud
ref_covm <- 0.01*diag(2)  # covariance matrix of the reference point-cloud  

dgts <- 6

# Generating the reference point-cloud

refOk <- F; set.seed(42) # this seed is used only for the generation of P_, which will be the same over all the independent runs!
while( !refOk ) {
  P_ <- rmvnorm( m, ref_mean, ref_covm )
  aux <- as.matrix(dist(P_))
  diag(aux) <- max(aux)
  refOk <- (min(aux)<=rho)
}


#****************************************************************************
# Random search in the Physical Space with N = N*nRuns of BO experiments
#****************************************************************************

bestSeens <- NULL
for( seed in 1:nSeeds ) {
  cat("> RS in Physcal Space:\n")
  set.seed(seed)
  vPs <- NULL; ys <- NULL
  for( i in 1:N ) {
    P <- randomLHS(m,2)
    vP <- as.vector(P)
    vPs <- rbind( vPs, vP)
    # check for feasibilty
    if( min(dist(P)>=rho) ) {
      tmp <- NULL
      for( j in 1:nfereps )
        # tmp <- c(tmp,objective1(vP))
        tmp <- c(tmp,permutation_invariant_objective1(vP))
      ys <- c( ys, min(tmp) )
    } else {
      ys <- c( ys, NA )
    }
    cat(".")
    if( i %%100 == 0 )
      cat(i,"/",N,"\n")
  }
  
  bestSeens <- c(bestSeens,min(ys,na.rm=T))
  
  cat("> Feasible solutions:",sum(!is.na(ys)),"/",N,"\n")
  cat("> Best y =",min(ys,na.rm=T),"\n")
  cat("> Repeated solutions in LHS ( with decimal digits =",dgts,"):",N-nrow(unique(data.frame(round(vPs,dgts)))),"\n")
  
  for( i in 1:N ) {
    P <- matrix(vPs[i,],m,2)
    ot <- transport( pp(P_), pp(P), p=2, method="networkflow" )
    stopifnot( nrow(ot)==m )
    ot <- ot[order(ot$from),]
    vPs[i,] <- as.vector(P[ot$to,])
  }
  cat("> Replica (i.e. permutations) w.r.t. to OT (with decimal digits =",dgts,"):",N-nrow(unique(data.frame(round(vPs,dgts)))),"\n")
  
}

bestSeens_physical <- bestSeens



#****************************************************************************
# Random search in the Flow Space with N = N*nRuns of BO experiments
#****************************************************************************


bestSeens <- NULL
for( seed in 1:nSeeds ) {
  cat("> RS in the Flow Space (i.e., GP's input space):\n")
  set.seed(seed)
  ys <- NULL
  vXs <- randomLHS(N,2*m)
  vPs <- NULL
  for( i in 1:N ) {
    vX <- as.vector(vXs[i,]-P_) #  rescale on the fly
    P <- P_+matrix(vX,m,2)
    vPs <- rbind( vPs, as.vector(P) )
    # check for feasibilty
    if( min(dist(P)>=rho) ) {
      tmp <- NULL
      for( j in 1:nfereps )
        # tmp <- c(tmp,objective1(as.vector(P)))
        tmp <- c(tmp,permutation_invariant_objective1(as.vector(P)))
      ys <- c( ys, min(tmp) )
    } else {
      ys <- c( ys, NA )
    }
    cat(".")
    if( i %%100 == 0 )
      cat(i,"/",N,"\n")
  }
  
  bestSeens <- c(bestSeens,min(ys,na.rm=T))
  
  cat("> Feasible solutions:",sum(!is.na(ys)),"/",N,"\n")
  cat("> Best y =",min(ys,na.rm=T),"\n")
  cat("> Repeated solutions in LHS ( with decimal digits =",dgts,"):",N-nrow(unique(data.frame(round(vXs,dgts)))),"\n")
  
  for( i in 1:N ) {
    P <- matrix(vPs[i,],m,2)
    ot <- transport( pp(P_), pp(P), p=2, method="networkflow" )
    stopifnot( nrow(ot)==m )
    ot <- ot[order(ot$from),]
    vPs[i,] <- as.vector(P[ot$to,])
  }
  cat("> Replica (i.e. permutations) w.r.t. to OT (with decimal digits =",dgts,"):",N-nrow(unique(data.frame(round(vPs,dgts)))),"\n")
  
  
}


cat("\n\n RS in physical space:",mean(bestSeens_physical),"[",min(bestSeens_physical),",",max(bestSeens_physical),"]\n")
cat("\n\n RS in flow space:",mean(bestSeens),"[",min(bestSeens),",",max(bestSeens),"]\n")
rm(list=ls()); graphics.off(); cat("\014")

cat("> Importing R packages..."); library(transport); library(DiceKriging); library(lhs); library(mvtnorm); library(reticulate); cat("done!\n")
cat("> Importing Python functions..."); source_python("ensemble_model_predict.py"); cat("done!\n")
cat("> Importing utils for permutation-invariant evaluation of the ensemble model..."); source("utils.R"); cat("done!\n")


#*********************************************************************************
# Experimental setting
#*********************************************************************************

# ***************************************************
# Problem setup
# ***************************************************

m <- 5        # number of turbines
rho <- 0.1512 # from 2 * rotor diam / farm length = (2*126)/(333.33*5)



# ***************************************************
# Bayesian Optimization setup
# ***************************************************

seeds <- 1:30           # seeds for independent runs
kernel <- "exp"         # GP's kernel (i.e., 'exp' and 'gauss' in the paper)
lcb.beta <- 1           # beta for GP's LCB (i.e., we minimize the amount of generated power changed in sign)
RSsamples <- 10000      # LHS samples for inexact optimization of the acquisition function (GP-LCB)  
n0 <- 2*m+1             # initial random solutions (i.e., 2*m+1 is the minimum required) to fit a GP)
N <- 500                # total number of queries (including initial random solutions)
nfereps <- 1            # repetitions of function evaluation for each query (the min value is taken)

ref_mean <- c(-1,-1)      # bi-variate mean of the reference point-cloud
ref_covm <- 0.01*diag(2)  # covariance matrix of the reference point-cloud  




#*********************************************************************************
# Main
#*********************************************************************************

# Generating the reference point-cloud
refOk <- F; set.seed(42) # this seed is used only for the generation of P_, which will be the same over all the independent runs!
while( !refOk ) {
  P_ <- rmvnorm( m, ref_mean, ref_covm )
  aux <- as.matrix(dist(P_))
  diag(aux) <- max(aux)
  refOk <- (min(aux)<=rho)
}
DP_ <- as.matrix(dist(P_)) # it will be used later

RES <-  data.frame( seed=rep(NA,N*length(seeds)),
                    iter=rep(NA,N*length(seeds)),
                    nFeasCandidates=rep(NA,N*length(seeds)),
                    x=matrix(NA,N*length(seeds),2*m),
                    acqValue_=rep(NA,N*length(seeds)),
                    y=rep(NA,N*length(seeds)),
                    trnTime=rep(NA,N*length(seeds)),
                    acqTime=rep(NA,N*length(seeds)),
                    evalTime=rep(NA,N*length(seeds)),
                    stringsAsFactors=F ) # data frame for storing results

rowIx <- 1

for( seed in seeds ) {
  
  cat("\n[ ***** SEED =",seed,"***** ]\n")
  set.seed(seed)
  
  # BO initialization
  cat("> BO initialization...")
  Xs <- NULL; ys <- numeric()
  for( i in 1:n0 ) {
    
    P <- maximinLHS(m,2) # LHS in the physical space (i.e., point-clouds)
    
    evalTime <- Sys.time()
    tmp <- numeric()
    for( k in 1:nfereps )
      tmp[k] <- permutation_invariant_objective(as.vector(P))
    evalTime <- difftime(Sys.time(),evalTime,units="secs")
    ys[i] <- min(tmp)
    
    # *****************************************************************
    # For vanilla BO on flows we do not need to retrieve the optimal
    # flows in terms of OT
    # *****************************************************************
    
    Xs <- rbind( Xs, as.vector(P-P_) )
    
    # storing into the data frame
    RES[rowIx,] <- data.frame( seed=seed,
                               iter=0,
                               nFeasCandidates=NA,
                               x=t(Xs[nrow(Xs),]),
                               acqValue_=NA,
                               y=ys[length(ys)],
                               trnTime=NA, acqTime=NA, evalTime=evalTime,
                               stringsAsFactors=F )
  }
  cat("done!\n")
  
  
  # BO sequential queries
  cat("> BO sequential query:\n  [")
  while( length(ys)<N ) {
    
    # GP fit
    trnTime <- Sys.time()
    gp <- km( design=data.frame(x=Xs), response=ys, covtype=kernel, nugget.estim=F, control=list(trace=0) )
    trnTime <- difftime(Sys.time(),trnTime,units="secs")
    
    # acquisition
    acqTime <- Sys.time()
    toSample <- T
    while( toSample ) {
      
      XX <- randomLHS(RSsamples,m*2)
      
      # rescaling XX w.r.t. to the original search space of P!
      for( j in 1:ncol(XX) )
        XX[,j] <- XX[,j]-as.vector(P_[j])
      
      # selecting only feasible X
      feasibleIxs <- numeric() 
      for( i in 1:nrow(XX) ) {
        X <- matrix(XX[i,],m,2)
        DX <- as.matrix(dist(X))
        ixs <- which(DX<DP_+rho, arr.ind=T)
        if( nrow(ixs)==m )
          feasibleIxs <- c(feasibleIxs,i)
      }
      toSample <- length(feasibleIxs)==0
      if( !toSample && length(feasibleIxs)<nrow(XX) ) {
        XX <- XX[feasibleIxs,]
      }
    }
    
    # optimizing (inexactly) the acquisition function
    aux <- predict( gp, data.frame(x=XX), "UK" )
    X_next <- XX[which.min(aux$mean-lcb.beta*aux$sd),]
    
    # making X_next consistent with OT!
    P <- P_ + matrix(X_next,m,2)

    
    # *****************************************************************
    # For vanilla BO on flows we do not need to retrieve the optimal
    # flows in terms of OT
    # *****************************************************************
    
    acqTime <- difftime( Sys.time(), acqTime, units="secs" )
    
    # function evaluation and update
    evalTime <- Sys.time()
    tmp <- numeric()
    for( k in 1:nfereps )
      tmp[k] <- permutation_invariant_objective(as.vector(P))
    evalTime <- difftime( Sys.time(), evalTime )
    
    Xs <- rbind( Xs, as.vector(P-P_) )
    ys <- c(ys,min(tmp))
    
    
    # storing into the data frame
    RES[rowIx,] <- data.frame( seed=seed,
                               iter=length(ys)-n0,
                               nFeasCandidates=nrow(XX),
                               x=t(Xs[nrow(Xs),]),
                               acqValue_=min(aux$mean-lcb.beta*aux$sd),
                               y=ys[length(ys)],
                               trnTime=trnTime, acqTime=acqTime, evalTime=evalTime,
                               stringsAsFactors=F )
    
    cat("=")
    if( (length(ys)-n0)%%50==0 || length(ys)==N ) {
      cat("] ",length(ys),"/",N,"\t(y+ = ",round(min(ys),2),")\n",sep="")
      if( length(ys)<N )
        cat("  [")
    }
  }
  cat("> best seen:",round(min(ys),6),"\n")
  
}

cat("\n\n> Summary of results:\n")
print(summary(aggregate(RES$y,by=list(RES$seed),min)$x))


BS <- NULL
for( s in sort(unique(RES$seed)) )
  BS <- rbind( BS, cummin( c(min(RES$y[RES$seed==s][1:n0]), RES$y[RES$seed==s][-(1:n0)]) ) )

BS_avg <- apply(BS,2,mean)
BS_up <- BS_avg + apply(BS,2,sd); BS_lo <- 2*BS_avg - BS_up

plot( 0:(N-n0), BS_avg, type="l", lwd=2, ylim=range(BS_lo, BS_up),
      xlab="BO iters", ylab="Best Seen" )
lines( 0:(N-n0), apply(BS,2,median), col="darkgrey", lwd=3 )
polygon( c(0:(N-n0),(N-n0):0),  c(BS_lo,rev(BS_up)), col=adjustcolor("black",alpha.f=0.1), border=F )

clrs <- rainbow(nrow(BS))
plot( 0:(N-n0), BS[1,], type="l", lwd=2, ylim=range(BS), col=clrs[1],
      xlab="BO iters", ylab="Best Seen" ) 
for( i in 2:nrow(BS) )
  lines( 0:(N-n0), BS[i,], type="l", lwd=2, col=clrs[i] ) 


cat("> Saving results...")
if( !dir.exists("RESULTS_vanilla_BO_on_flows") )
  dir.create("RESULTS_vanilla_BO_on_pointclouds")
saveRDS( RES, paste0("RESULTS_vanilla_BO_on_flows/results_",length(seeds),"_",kernel,"_",lcb.beta,"_",n0,"_",N,"_",nfereps,".RDS") )
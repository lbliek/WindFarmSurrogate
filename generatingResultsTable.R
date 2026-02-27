rm(list=ls()); graphics.off(); cat("\014")

folders <- c("RESULTS", "RESULTS_vanilla","RESULTS_vanilla_physical")

for( folder in folders ) {
  
  if( folder=="RESULTS" ) {
    cat(" [ PIDO ]\n")
  } else {
    if( folder=="RESULTS_vanilla" ) {
      cat(" [vannila BO on Flows]\n")
    } else {
      cat(" [vanilla BO on point-clouds (Physical Space) ]\n")
    }
  }
  
  allFiles <- list.files(folder)
  
  TABLE <- NULL
  for( f in allFiles ) {
    res <- readRDS( paste0(folder,"/",f) )
    
    aggr <- aggregate(res$y,by=list(res$seed),min)
    
    info <- unlist(strsplit(f,"_",fixed=F))
    nSeeds <- as.numeric(info[2])
    kernel <- info[3]
    lcb.beta <- as.numeric(info[4])
    n0 <- as.numeric(info[5])
    N <- as.numeric(info[6])
    nfeReps <- as.numeric(info[7])
    pushTheLimit <- info[8]=="TRUE.RDS"
    
    TABLE <- rbind( TABLE, data.frame( nSeeds=nSeeds,
                                       kernel=kernel,
                                       lcb.beta=lcb.beta,
                                       n0=n0,
                                       N=N,
                                       nfeReps=nfeReps,
                                       pushTheLimit=pushTheLimit,
                                       bs_med=round(median(aggr$x),2),
                                       bs_avg=round(mean(aggr$x),2),
                                       bs_min=round(min(aggr$x),2),
                                       bs_max=round(max(aggr$x),2),
                                       stringsAsFactors=F ))
  }
  
  print(TABLE)
  cat("\n\n")
}



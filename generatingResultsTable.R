rm(list=ls()); graphics.off(); cat("\014")

folders <- c("RESULTS_PIBO", "RESULTS_vanilla_BO_on_flows","RESULTS_vanilla_BO_on_pointclouds")

for( folder in folders ) {
  
  if( folder=="RESULTS_PIBO" ) {
    cat(" [ PIDO ]\n")
  } else {
    if( folder=="RESULTS_vanilla_BO_on_flows" ) {
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
    
    TABLE <- rbind( TABLE, data.frame( nSeeds=nSeeds,
                                       kernel=kernel,
                                       lcb.beta=lcb.beta,
                                       n0=n0,
                                       N=N,
                                       bs_avg=round(mean(aggr$x),2),
                                       bs_sd=round(sd(aggr$x),2),
                                       bs_med=round(median(aggr$x),2),
                                       bs_min=round(min(aggr$x),2),
                                       bs_max=round(max(aggr$x),2),
                                       stringsAsFactors=F ))
  }
  
  print(TABLE)
  cat("\n\n")
}






for( f in allFiles ) {
  r1 <- readRDS(paste0(folders[1],"/",f) )
  r2 <- readRDS(paste0(folders[2],"/",f) )
  r3 <- readRDS(paste0(folders[3],"/",f) )
  r1 <- aggregate(r1$y,by=list(r1$seed),min)
  r2 <- aggregate(r2$y,by=list(r2$seed),min)
  r3 <- aggregate(r3$y,by=list(r3$seed),min)
  aux <- cbind(r1$x,r2$x,r3$x)
  aux <- apply(aux,1,which.min)
  cat(f,"  \t won by PIBO:",sum(aux==1),"/",nrow(r1),"\t won by BO-on-flows:",sum(aux==2),"/",nrow(r1),"\n")
}
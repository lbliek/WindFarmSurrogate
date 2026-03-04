ext <- read.table("TPE_RS_WWsorted_surrogate_iters.csv",header=T,sep="," )

tpe <- ext[ext$approach=="hyperopt/tpe",]

ys <- numeric(nrow(tpe))
for( i in 1:nrow(tpe) ) {
  aux <- gsub(pattern="np.float64(",replacement="",tpe$iter_x[i],fixed=T)
  aux <- gsub(")","",aux,fixed=T); aux <- gsub("[","",aux,fixed=T); aux <- gsub("]","",aux,fixed=T)
  aux <- as.numeric(unlist(strsplit(aux,", ",fixed=T)))
  ys[i] <- permutation_invariant_objective(aux)
}
tpe <- cbind(tpe,ys)

tpebs <- NULL
for( exp_id in sort(unique(tpe$exp_id)) ) {
  aux <- tpe$ys[tpe$exp_id==exp_id]
  tpebs <- rbind(tpebs,cummin(c(min(aux[1:(2*m+1)]),aux[-(1:(2*m+1))])))
}
tpem <- -apply(tpebs,2,mean); tpelo <- tpem - apply(tpebs,2,sd); tpeup <- tpem + apply(tpebs,2,sd)

polygon( c(0:(N-n0),(N-n0):0), c(tpelo,rev(tpeup)), col=adjustcolor("red",alpha.f=0.3), border=F )
lines( 0:(N-n0), tpem, col="red", lwd=3 )

legend("bottomright",legend=c("PIBO","TPE"), col=c("gray","red"), lwd=4 )


stop(".....")

tpe <- ext[ext$approach=="hyperopt/tpe",-c(8,10)]


tpebs <- NULL
for( exp_id in sort(unique(tpe$exp_id)) ) {
  aux <- tpe$iter_best_fitness[tpe$exp_id==exp_id]
  tpebs <- rbind(tpebs,aux[-(1:(2*m))])
}
tpem <- -apply(tpebs,2,mean); tpelo <- tpem - apply(tpebs,2,sd); tpeup <- tpem + apply(tpebs,2,sd)

polygon( c(0:(N-n0),(N-n0):0), c(tpelo,rev(tpeup)), col=adjustcolor("red",alpha.f=0.3), border=F )
lines( 0:(N-n0), tpem, col="red", lwd=3 )

legend("bottomright",legend=c("PIBO","TPE"), col=c("gray","red"), lwd=4 )

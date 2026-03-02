ext <- read.table("TPE_RS_WW_sorted_surrogate_iters.csv",header=T,sep="," )

tpe <- ext[ext$approach=="hyperopt/tpe",-c(8,10)]

tpebs <- NULL
for( exp_id in sort(unique(tpe$exp_id)) ) {
  aux <- tpe$iter_best_fitness[tpe$exp_id==exp_id]
  tpebs <- rbind(tpebs,aux[-(1:(2*m))])
}
tpem <- -apply(tpebs,2,mean); tpelo <- tpem - apply(tpebs,2,sd); tpeup <- tpem + apply(tpebs,2,sd)

polygon( c(0:(N-n0),(N-n0):0), c(tpelo,rev(tpeup)), col=adjustcolor("pink",alpha.f=0.3), border=F )
lines( 0:(N-n0), tpem, col="red", lwd=3 )

legend("bottomright",legend=c("PIBO","TPE"), col=c("gray","red"), lwd=4 )

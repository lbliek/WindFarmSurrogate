rm(list=ls()); graphics.off(); cat("\014")
source("utils.R")
library(reticulate); source_python("ensemble_model_predict.py")
m <- 5
N <- 500
n0 <- 2*m+1
RES <- readRDS("RESULTS/results_30_exp_1_11_500_1.RDS")
bs <- NULL
for( i in sort(unique(RES$seed)) ) {
  tmp <- RES[RES$seed==i,]
  bs <- rbind( bs, cummin(c(min(tmp$y[tmp$iter==0]),tmp$y[tmp$iter!=0])) )
}
bsm <- -apply(bs,2,mean); bsl <- bsm-apply(bs,2,sd); bsu <- bsm+apply(bs,2,sd)
plot( NA, NA, xlim=c(0,N-n0), ylim=range(bsu,bsl), xlab="iters", ylab="f(x)"  )
polygon( c(0:(N-n0),(N-n0):0), c(bsl,rev(bsu)), col=adjustcolor("deepskyblue",alpha.f=0.2), border=F )
lines( 1:length(bsm)-1, bsm, lwd=3, col="blue" )




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

legend("bottomright",legend=c("PIBO","TPE"), col=c("blue","red"), lwd=4 )


# stop(".....")
# 
# tpe <- ext[ext$approach=="hyperopt/tpe",-c(8,10)]
# 
# 
# tpebs <- NULL
# for( exp_id in sort(unique(tpe$exp_id)) ) {
#   aux <- tpe$iter_best_fitness[tpe$exp_id==exp_id]
#   tpebs <- rbind(tpebs,aux[-(1:(2*m))])
# }
# tpem <- -apply(tpebs,2,mean); tpelo <- tpem - apply(tpebs,2,sd); tpeup <- tpem + apply(tpebs,2,sd)
# 
# polygon( c(0:(N-n0),(N-n0):0), c(tpelo,rev(tpeup)), col=adjustcolor("red",alpha.f=0.3), border=F )
# lines( 0:(N-n0), tpem, col="red", lwd=3 )
# 
# legend("bottomright",legend=c("PIBO","TPE"), col=c("gray","red"), lwd=4 )

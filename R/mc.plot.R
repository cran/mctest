mc.plot<-function(x, y, Inter = FALSE, vif = 10, ev = 0.01, ...){

  #from lm.fit (extra argument handling)
  if(length(list(...))>1L)
  {warning("Extra arguments ", paste(sQuote(names(list(...) ) ) , sep=", "),
           " are ignored", domain=NA)}
  else if (length(list(...))==1L) warning("Extra argument ", sQuote(names(list(...) ) ),
                                          " is ignored", domain=NA)

  par(mfrow=c(2,1),mar=c(3,4,3,1))

  xcol=colnames(x)

  if(Inter==TRUE || is.null(Inter) ){
    Eigval<-eigprop(x, Inter = TRUE)$ev
    xs=length(xcol)+1
    lends=c("Intercept", xcol)
  }else{
    Eigval<-eigprop(x, Inter = FALSE)$ev
    xs=length(xcol)
    lends=xcol
  }

  VIF<-imcdiag(x, y, method = "VIF")[[1]][,1]
# VIF Plot and setting
  plot(VIF,
       lty=1,
       type="h",
       xaxt="n",
       main="VIF Plot",
       ylab="VIF values",
       xlab="",
       ylim=c(0, max(VIF)),
       lwd=2
      )
  axis(side=1,
       at=1:length(xcol),
       labels=colnames(x),
       las=0
       )
  text(1:length(xcol),
       VIF,
       round(VIF,3),
       col="blue",
       cex=.95
       )
  abline(h=vif,
         col="red",
         lty=2,
         lwd=2
         )
  text(1,
       vif,
       paste("vif=",vif),
       col="blue",
       cex=1,
       pos=4
       )
#  mtext("Regressors",  side=1,padj=-2, outer=TRUE)
  #title("Multcollinearity Detection Plots")
  #points(vif,cex=1, col="red")

# Eigenvalues plot and settings
  plot(Eigval,
       lty=1,
       type="l",
       col=1:length(lends),
       main="Eigenvalues Plot",
       xaxt="n",
       ylab="Eigenvalues",
       xlab="",
       lwd=2
       )
  axis(side=1,
       at=1:xs,
       labels=paste("EV", 1:xs),
       las=0
       )
  text(1:xs,
       Eigval,
       round(Eigval,3),
       col="blue",
       cex=.95
       )
  abline(h=ev,
         col="red",
         lty=2,
         lwd=2
         )
  text(1,
       ev,
       paste("EV=",ev),
       col="blue",
       cex=1,
       pos=4
       )

  #mtext("Regressors",  side=1,padj=-2, outer=TRUE)
}

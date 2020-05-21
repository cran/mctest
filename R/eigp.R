eigprop <- function(mod, na.rm = TRUE, Inter = TRUE, prop = 0.5, ...){

  cl<-match.call()
  if (!is.null(mod$call$formula)){
    if (Inter==FALSE)
    x <- as.matrix(model.frame(mod)[,-1]) # Regressors only 
    #y <- as.matrix(model.frame(mod)[,1]) # dependent variable
  }
  if (Inter==TRUE) 
    x <- as.matrix(mod$model)
  
#  x<-as.matrix(x)
#  y<-as.matrix(y)

  match.call()

  #from lm.fit (extra argument handling)
  if(length(list(...))>1L)
  {
    warning("Extra arguments ", paste(sQuote(names(list(...) ) ) , sep=", "),
           " are ignored", domain=NA)
  }
  else if (length(list(...))==1L) warning("Extra argument ", sQuote(names(list(...) ) ),
                                          " is ignored", domain=NA)

  if(ncol(x) < 2)
    stop('X matrix must contain more than one variable')

#  if(!is.numeric(x) || !is.numeric(y))
  if(!is.numeric(x))
    stop('X must be a numeric matrix')

#  if(nrow(x)!=length(y))
#    stop('X and y contain different numbers of observations')

  #remove the missing values and re-create the data set
  if( na.rm ) {
    df<-as.data.frame(x) #data
#    ncolxy<-ncol(df)
    df<-df[complete.cases(df), ]
#  x<-na.omit(x)
#    x<-as.matrix(x)
#    y<-as.matrix(df[,ncolxy])
    x<-as.matrix(df)
 }

#  nvar<-ncol(x)
#  n<-nrow(x)

#if(Inter==TRUE){
#  ex<-cbind(1,x)
#  colnames(ex)<-c("Intercept", colnames(x))
#}else {
#  ex<-x
#  colnames(ex)<-colnames(x)
#}

xz<-apply(x,2,function(x){x/sqrt(sum(x^2))})

corxz<-t(xz)%*%xz

#ev<-.Internal(La_rs(corxz, FALSE))$values;ev
ev<-eigen(corxz)$values
#sdev<-sqrt(ev)
#prop<-ev/sum(ev)
#cum<-cumsum(ev)/sum(ev)
svdX<-svd(xz)
ordev<-sort(ev, decreasing = T)
ci<-sqrt(max(ev)/ev)
phi=t((svdX$v%*%diag(1/sqrt(ordev)))^2)
#phi=svdX$v%*%diag(1/svdX$d);phi
#phi=t(phi^2);phi
pi<-prop.table(phi,2)
colnames(pi)<-colnames(xz)
rownames(pi)<-1:nrow(pi)

eigpres<-list(ev = ev,
              ci = ci,
              pi = pi,
              call = cl,
              Inter = Inter,
              prop = prop)
class(eigpres) <- "eigp"
eigpres
}

print.eigp<-function(x, digits = max(3, getOption("digits") - 3), ...){

  res<-cbind(x$ev,
             x$ci,
             x$pi)

  colnames(res)<-c("Eigenvalues",
                   "CI",
                   colnames(x$pi))

  rownames(res)<-1:nrow(res)

  fcall=cat("\nCall:\n",
            paste(deparse(x$call),
                  sep = "\n",
                  collapse = "\n"),
            "\n\n", sep = "")

  fcall

  print(round(res, digits))

  if(x$prop >= 1.0){
    cat("\n==================Note===============\n")
    cat("Variance Proportion threshold used is greater than 1\n")
    }
  else {
  if(x$Inter==TRUE){
    pi<-x$pi[,-1]
    ix <- pi >= x$prop
    if(all(ix!=1)){
      cat("\n===============================")
      cat(sprintf("\nnone of the variance proportion is > %0.02f", x$prop),"\n" )
    }
      else{
      cat("\n===============================")
      cat(sprintf("\nRow %s==> %s, proportion %f >= %0.02f",
              rownames(pi)[row(pi)[ix]],
              colnames(pi)[col(pi)[ix]],
              pi[ix], x$prop),"\n")
      }
  }
  else {
    pi<-x$pi
    ix<- pi >= x$prop
    if(all(ix!=1)){
      cat("\n===============================")
      cat(sprintf("\nnone of the variance proportion is > %0.02f", x$prop), "\n")
    }
    else{
    cat("\n===============================")
    cat(sprintf("\nRow %s==> %s, proportion %f >= %0.02f",
                rownames(x$pi)[row(x$pi)[ix]],
                colnames(x$pi)[col(x$pi)[ix]],
                x$pi[ix], x$prop ),"\n")
    }
  }
  }

}

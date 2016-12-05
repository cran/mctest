omcdiag<-function(x, y, na.rm=TRUE, Inter=TRUE, detr=0.01, red=0.5,
                  conf=0.95, theil=0.5,cn=30,...){
  cl<-match.call()

  x<-as.matrix(x)
  y<-as.matrix(y)

  match.call()

  #from lm.fit (extra argument handling)
  if(length(list(...))>1L)
  {warning("Extra arguments ", paste(sQuote(names(list(...) ) ) , sep=", "),
           " are ignored", domain=NA)}
  else if (length(list(...))==1L) warning("Extra argument ", sQuote(names(list(...) ) ),
                                          " is ignored", domain=NA)

  if(ncol(x)<2)
    stop('X matrix must contain more than one variable')

  if(!is.numeric(x) || !is.numeric(y))
    stop('X must be a numeric matrix')

  if(nrow(x)!=length(y))
    stop('X and y contain different numbers of observations')

  #remove the missing values and re-create the data set
  if( na.rm ) {
    df<-as.data.frame(cbind(x,y)) #data
    ncolxy<-ncol(df)
    df<-df[complete.cases(df),]
    y<-as.matrix(df[,ncolxy])
    x<-as.matrix(df[,-ncolxy])
  }

  nvar<-ncol(x)
  n<-nrow(x)

  R2<-matrix(NA, nvar,1)
  AdjR2<-matrix(NA, nvar,1)
  R2yonx<-matrix(NA, nvar,1)

  #overall R-Square from all regressors
  R2yonallx<-summary(lm(y~x))$r.squared

  #R-square, adjusted-R-square and R-square due to single x (all from auxiliary regression)
  for(i in 1:nvar)   {
    R2[i]<-summary(lm(x[,i]~x[,-i]))$r.squared  #R-Square from Auxuliary Regression
    AdjR2[i]<-summary(lm(x[,i]~x[,-i]))$adj.r.squared #Adjusted R=square from Auxuliary Regression

    R2yonx[i]<-summary(lm(y~x[,i]))$r.squared
  }

  sx<-scale(x)/sqrt(n-1)
  R<-t(sx)%*%sx      # Correlation matrix from Design matrix
  Det<-det(R)        # Determinant of correlation matrix
  Det1<-cbind(Det, Det<detr)

  Fchi<- -(nrow(x)-1-(1/nvar) * (2*nvar + 5) ) * log(Det)
  Fchi<-cbind(Fchi, Fchi>qchisq(conf,1/2*(nvar)*(nvar-1)))

  Red <- sqrt((sum((t(sx)%*%sx)^2)-nvar)/(nvar*(nvar-1)))
  Red <-cbind(Red, Red>red)
  #sum of lambda^-1 values
  slambda<-sum(1/eigen(R)$values)
  slambda<- cbind(slambda, slambda>5*nvar)

  #Theil indicator
  Theil<-R2yonallx-sum(R2yonallx-R2)
  Theil<-cbind(Theil, Theil>theil)
  #R2yonallx<-summary(lm(y~x))$r.squared

  if(Inter==T){
    ex<-cbind(1,x)
  }else {
    ex<-x
    colnames(ex)<-names(x)
  }

  xz<-apply(ex,2,function(ex){ex/sqrt(sum(ex^2))})
  Eigval<-eigen(t(xz)%*%xz)$values
#  P<- eigen(t(xz)%*%xz)$vectors
  CN<-sqrt(max(Eigval)/min(Eigval))
  CN<-cbind(CN,CN>cn)

  odiags<-list(Det=Det1, Fchi=Fchi, Red=Red, slambda=slambda, Theil=Theil, CN=CN)
  odiags<-  do.call(rbind,odiags)
  colnames(odiags)<-c("results", "detection")
  rownames(odiags)<-c("Determinant", "Farrar Chi-Square", "Red Indicator",
                      "sum of Lambda Invers", "Theil Indicator", "Condition Number")
  ores<-list(odiags=odiags,
           #  nvar=nvar,
             Eigval=Eigval,
             Inter=Inter,
             x=x, calll=cl)

  ores<-c(ores)
  class(ores)<-"omc"
  ores
}

print.omc<-function(x,digits=max(3,getOption("digits")-3),...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("\nOverall Multicollinearity Diagnostics\n\n")
#  nvar=ncol(x$x)

  EigCI<-rbind(x$Eigval, sqrt(max(x$Eigval)/x$Eigval))
  rownames(EigCI)<-c("Eigenvalues:", "Condition Indeces:")

  omc<-round(x$odiags[,1],digits)
  omd<-x$odiags[,2]
  res<-cbind(omc, omd)
  rownames(res)<-c("Determinant |X'X|:",
                   "Farrar Chi-Square:",
                   "Red Indicator:",
                   "Sum of Lambda Inverse:",
                   "Theil's Method:",
                   "Condition Number:")
  colnames(res)<-c("MC Results", "detection")

  print(res)
  cat("\n1 --> COLLINEARITY is detected \n0 --> COLLINEARITY in not detected by the test\n\n")
  cat("===================================\n")

  if(x$Inter){
    cat("Eigvenvalues with INTERCEPT\n")
    colnames(EigCI)<-c("Intercept", colnames(x$x))
    print(round(EigCI,digits))
  }else{
    cat("Eigenvalues without INTERCEPT\n\n")
    colnames(EigCI)<-colnames(x$x)
    print(round(EigCI, digits))
  }
  invisible(res)
}

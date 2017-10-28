imcdiag<-function(x, y, method=NULL, na.rm=TRUE, corr=FALSE, vif=10, tol=0.1,
                  conf=0.95,cvif=10, leamer=0.1, all=FALSE,...){

  x<-as.matrix(x)
  y<-as.matrix(y)

  cl<-match.call()

  #from lm.fit (extra argument handling)
  if(length(list(...))>1L)
  {warning("Extra arguments ", paste(sQuote(names(list(...) ) ) , sep=", "),
           " are ignored", domain=NA)}
  else if (length(list(...))==1L) warning("Extra argument ", sQuote(names(list(...) ) ),
                                          " is ignored", domain=NA)

  if(is.null(corr) || corr==FALSE){
    corr=FALSE
  } else corr=TRUE


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

  nvar<-ncol(x)# Variables in Design matrix
  n<-nrow(x)     # Observations in data
  sx<-scale(x)/sqrt(n-1)

  R2<-matrix(NA, nvar,1)
  AdjR2<-matrix(NA, nvar,1)
  R2yonx<-matrix(NA, nvar,1)

  R2yonallx<-summary(lm(y~x))$r.squared

  #R-square, adjusted-R-square and R-square due single x (all from auxiliary regression)
  for(i in 1:nvar)   {
    R2[i]<-summary(lm(x[,i]~x[,-i]))$r.squared  #R-Square from Auxuliary Regression
    AdjR2[i]<-summary(lm(x[,i]~x[,-i]))$adj.r.squared #Adjusted R=square from Auxuliary Regression

    R2yonx[i]<-summary(lm(y~x[,i]))$r.squared
  }

  #pvalues of coefficients from OLS
  pval<-which(summary(lm(y~x))$coefficients[-1,4]>1-conf)

  #VIF
  VIF<-1/(1-R2)

  #TOL
  TOL<-1/VIF

  #Farrar wi test of F
  Wi<- (R2/(1-R2))*( (n-nvar)/(nvar-1) )

  #F and R-square relaion
  #Auxiliary F test Gujarati page 361
  #Fi<-(R2/(nvar-2))/((1-R2)/(n-nvar+1))
  Fi<-(R2/(nvar-2))/((1-R2)/(n-nvar+1))

  #leamer's Method
  matxdev<-((apply(sx,2,sd))^2)*(nrow(sx)-1)
  Leamer<-sqrt( (1/matxdev)/VIF) # Leamer's Calculations

  #CVIF
  CVIF<-VIF* (1-R2yonallx)/(1-sum(R2yonx))

  idiags<-vector("list")

  alldiag<-cbind(VIF>=vif,
                TOL<=tol,
                Wi>=qf(conf, n-nvar, nvar-1),
                Fi>=qf(conf, n-2, n-nvar+1),
                Leamer<=leamer,
                CVIF>=cvif,
                R2>R2yonallx
                )

  colnames(alldiag)<-cbind("VIF",
                           "TOL",
                           "Wi",
                           "Fi",
                           "Leamer",
                           "CVIF",
                           "Klein")
  rownames(alldiag)<-colnames(x)

  if(is.null(method)){
    idiags<-list(VIF=VIF,
                 TOL=TOL,
                 Wi=Wi,
                 Fi=Fi,
                 Leamer=Leamer,
                 CVIF=CVIF,
                 Klein=R2>R2yonallx)
    idiags<-do.call(cbind,idiags)
    colnames(idiags)<-cbind("VIF",
                            "TOL",
                            "Wi",
                            "Fi",
                            "Leamer",
                            "CVIF",
                            "Klein")
    rownames(idiags)<-colnames(x)
  }else  {
    idiags<-switch(method,
                VIF   =cbind(VIF, VIF>=vif),
                TOL   =cbind(TOL, TOL<=tol),
                Wi    =cbind(Wi,  Wi>=qf(conf, n-nvar, nvar-1)), # F(n-p, p-1)
                Fi    =cbind(Fi,  Fi>=qf(conf, n-2, n-nvar+1)),    # F(n-2, n-p+1)
                Klein =cbind(R2,  R2yonallx, R2-R2yonallx, R2>R2yonallx),
                Leamer=cbind(Leamer, Leamer<=leamer),
                CVIF  =cbind(CVIF, CVIF>=cvif),
                stop("\n\nThe argument of method should be VIF, TOL, Wi, CVIF, Klean, Leamer, or Fi\n\n")
    )
  }

  if (ncol(idiags)==2){
    colnames(idiags)<-c(method, "detection")
    rownames(idiags)<-colnames(x)
  } else if (ncol(idiags)==4){
    colnames(idiags)<-c("R2j",
                        "R2(overall)",
                        "Difference",
                        "detection")
    rownames(idiags)<-colnames(x)
  }

  ires<-list(idiags=idiags,
             x=x,
             y=y,
             method=method,
             corr=corr,
             call=cl,
             pval=pval,
             R2=R2yonallx,
             all=all,
             alldiag=alldiag)

  class(ires)<-"imc"
  ires
}

print.imc<-function(x, digits = max(3, getOption("digits") - 3), ...){
  n<-nrow(x$x)
  #nvar<-ncol(x$x)
  method<-x$method
  res<-x$idiags

  fcall=cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  if(is.null(method) && x$all==FALSE){
    fcall
    cat("\nAll Individual Multicollinearity Diagnostics Result\n\n")
    print(round(res,digits))}
  else if (!is.null(method) && x$all==FALSE) {
    fcall
    cat("\n", method, "Multicollinearity Diagnostics\n\n")
    print(round(res, digits))
  }else if (!is.null(method) && x$all==TRUE){
    fcall
    cat("\n", method, "Multicollinearity Diagnostics\n\n")
    print(round(res, digits ))
    cat("\nAll Individual Multicollinearity Diagnostics in 0 or 1 \n\n")
    print(1*x$alldiag)
  } else{
    fcall
    cat("\nAll Individual Multicollinearity Diagnostics in 0 or 1 \n\n")
    print(1*x$alldiag)
  }

   #cat("\n")

  if(!is.null(method) && ncol(res)==2 && sum(res[,2]!=0) ){
    cat("\nMulticollinearity may be due to", colnames(x$x) [which(res[,2] %in% TRUE)],"regressors\n")
    cat("\n1 --> COLLINEARITY is detected by the test \n0 --> COLLINEARITY is not detected by the test\n\n")
    cat("===================================\n")
  } else if (!is.null(method) && ncol(res)==4 && sum(res[,4]!=0)){
    cat("\nMulticollinearity may be due to", colnames(x$x) [which(res[,4] %in% TRUE)],"regressors\n\n")
    cat("\n1 --> COLLINEARITY is detected by the test \n0 --> COLLINEARITY is not detected by the test\n\n")
    cat("===================================\n")
  } else if (!is.null(method) ){
    cat("\nNOTE: ", method, "Method Failed to detect multicollinearity\n\n")
    cat("\n0 --> COLLINEARITY is not detected by the test\n\n")
    cat("===================================\n")
  } else {cat("\n1 --> COLLINEARITY is detected by the test \n0 --> COLLINEARITY is not detected by the test\n\n")
    #cat(x$pval, ", have non-significant t-ratios\n")
    if(!length(x$pval)){
      cat("* all coefficients have significant t-ratios\n")
    } else {
      cat(paste(colnames(x$x)[x$pval],","), "coefficient(s) are non-significant may be due to multicollinearity\n")
    } #colnames(x$x)[x$pval] )

    cat("\nR-square of y on all x:", round(x$R2, digits),"\n")
    cat("\n* use method argument to check which regressors may be the reason of collinearity\n")
    cat("===================================\n")
  }

  if(x$corr==TRUE){
    cat("\nCorrelation Matrix\n")
    sx<-scale(x$x)/sqrt(n-1)
    corR<-t(sx)%*%sx
    ix <- abs(corR) >= 0.7 & upper.tri(corR)

    print(corR)
    cat("\n====================NOTE===================\n")
    cat(sprintf("\n%s and %s may be collinear as |%f|>=0.7", colnames(corR)[row(corR)[ix]],
                colnames(corR)[col(corR)[ix] ], corR[ix] ),"\n\n")
  }

}

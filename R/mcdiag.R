mctest<-function(mod, type=c("o","i","b"), na.rm = TRUE, Inter=TRUE, method=NULL, corr=FALSE,
                 detr=0.01, red=0.5, theil=0.5, cn=30, vif=10, tol=0.1, conf=0.95, cvif=10,
                 ind1=0.02, ind2=0.7, leamer=0.1,all=FALSE,...){

  type<-match.arg(type)

    if (!is.null(mod$call$formula)){
      x <- model.matrix(mod)[, -1] # Regressors only 
      y <- as.matrix(model.frame(mod)[, 1]) # dependent variable
  }
  

    if(ncol(x)<2)
        stop('X matrix must contain more than one variable')

#    if(!is.numeric(x) || !is.numeric(y))
#        stop('X must be a numeric matrix')

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


    match.call()

    #from lm.fit (extra argument handling)
    if(length(list(...))>1L)
    {warning("Extra arguments ", paste(sQuote(names(list(...) ) ) , sep=", "),
             " are ignored", domain=NA)}
    else if (length(list(...))==1L) warning("Extra argument ", sQuote(names(list(...) ) ),
                                            " is ignored", domain=NA)


    if(type=="o" ){
      if(!missing(corr) || !missing(method))
        warning("\n\nmethod or corr argument is required for imcdiag function only\n")
      if (Inter==FALSE)
        (omcdiag(mod, Inter=FALSE,detr=detr, red=red, conf=conf, theil=theil,cn=cn,...))
      else (omcdiag(mod,Inter=TRUE, detr=detr, red=red, conf=conf, theil=theil,cn=cn,...))

    } else if (type=="b"){
      if(Inter==FALSE || corr==FALSE || is.null(corr) || is.null(Inter) || is.null(all) || all==FALSE){
        print(omcdiag(mod, Inter=Inter, detr=detr, red=red, conf=conf, theil=theil,cn=cn,...))
        cat("\n===================================\n")
        print(imcdiag(mod,method,corr=FALSE, vif=vif, tol=tol, conf=conf, cvif=cvif,
                      ind1=ind1, ind2=ind2, leamer=leamer, all=all))
      }

      else  {
        print(omcdiag(mod,Inter=Inter,detr=detr, red=red, conf=conf, theil=theil,cn=cn,...))
        cat("\n===================================\n")
        print(imcdiag(mod,method, corr=TRUE, vif=vif, tol=tol, conf=conf, cvif=cvif,
                      ind1=ind1, ind2=ind2, leamer=leamer, all=all))
      }

      }
    else{
      if(!missing(Inter))
        warning("\n\nInter argument is required for omcdiag function\n")
      if (corr==FALSE || is.null(corr))
        (imcdiag(mod, method,corr=FALSE,vif=vif,tol=tol, conf=conf,cvif=cvif,
                 ind1=ind1, ind2=ind2, leamer=leamer, all=all) )
      else (imcdiag(mod, method, corr=TRUE,vif=vif,tol=tol, conf=conf,cvif=cvif,
                    ind1=ind1, ind2=ind2, leamer=leamer, all=all) )
    }
}

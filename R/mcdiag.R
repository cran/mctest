mctest<-function(x, y, type=c("o","i","b"), na.rm = TRUE, Inter=TRUE, method=NULL, corr=FALSE,
                 detr=0.01, red=0.5, theil=0.5,cn=30,vif=10, tol=0.1, conf=0.95, cvif=10,
                 leamer=0.1,...){

    x<-as.matrix(x)
    y<-as.matrix(y)

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

#    n<-nrow(x)

    if(type=="o" || is.null(type)){
        if(!missing(corr) || !missing(method))
            warning("\n\nmethod or corr argument is required for imcdiag function only\n")
        if (Inter==FALSE)
             (omcdiag(x,y,Inter=FALSE,detr=detr, red=red, conf=conf, theil=theil,cn=cn,...))
        else (omcdiag(x,y,Inter=TRUE, detr=detr, red=red, conf=conf, theil=theil,cn=cn,...))

    }else if (type=="b"){
        if(Inter==FALSE || corr==FALSE || is.null(corr) || is.null(Inter)){
            print(omcdiag(x,y, Inter=Inter, detr=detr, red=red, conf=conf, theil=theil,cn=cn,...))
            cat("\n===================================\n")
            print(imcdiag(x,y,method,corr=FALSE, vif=vif,tol=tol, conf=conf,cvif=cvif,
                          leamer=leamer))
            }
        else {
            print(omcdiag(x,y,Inter=Inter,detr=detr, red=red, conf=conf, theil=theil,cn=cn,...))
            cat("\n===================================\n")
            print(imcdiag(x,y,method,corr=TRUE,vif=vif,tol=tol, conf=conf,cvif=cvif,
                          leamer=leamer))
        }
    }
    else{
        if(!missing(Inter))
            warning("\n\nInter argument is required for omcdiag function\n")
        if (corr==FALSE || is.null(corr))
           (imcdiag(x, y, method,corr=FALSE,vif=vif,tol=tol, conf=conf,cvif=cvif,
                    leamer=leamer) )
        else (imcdiag(x, y, method,corr=TRUE,vif=vif,tol=tol, conf=conf,cvif=cvif,
                    leamer=leamer) )
    }
}

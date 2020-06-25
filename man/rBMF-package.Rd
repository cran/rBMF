\name{rBMF-package}
\alias{rBMF-package}
\alias{rBMF}
\encoding{utf8}
\docType{package}
\title{
\packageTitle{rBMF}
}
\description{
\packageDescription{rBMF}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{rBMF}
\packageIndices{rBMF}

}
\author{Abdelmoneim Amer Desouki}

\references{
   
   topFiberM
   -Desouki, A. A., Röder, M., & Ngomo, A. C. N. (2019). 
        topFiberM: Scalable and Efficient Boolean Matrix Factorization.
        arXiv preprint arXiv:1903.10326.
   
   Asso
   -Miettinen, P., Mielikäinen, T., Gionis, A., Das, G., & Mannila, H. (2008). 
    The discrete basis problem. IEEE transactions on knowledge and data engineering,
    20(10), 1348-1362.
   
   GreConD, GreConDPlus
   -Belohlavek R., Vychodil V.: 
     Discovery of optimal factors in binary data via a novel method of matrix decomposition. 
     Journal of Computer and System Sciences 76(1)(2010), 3-20
}
\keyword{ Boolean Matrix Factorization }

\seealso{
%~~ Optional links to other man pages, e.g. ~~
%\code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
    \code{\link{topFiberM}} 
    \code{\link{Asso_approximate}} 
    \code{\link{GreConD}} 
    \code{\link{GreConDPlus}}
}
\examples{
%~~ simple examples of the most important functions ~~
\donttest{
data(DBLP)
 X=DBLP
    r=7
     Xb=X==1#Convert to boolean
    tempX=as(X,'TsparseMatrix')
    stats=NULL
    for(tP in c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,1)){
    
      Res=topFiberM(Xb,r=r,tP=tP,SR=100,verbose=1)
    
    X_=Res$A \%*\% Res$B
    X_=as(X_,'TsparseMatrix')
    #Calculate metrics
    li=tempX@i[tempX@x==1]+1
    lj=tempX@j[tempX@x==1]+1
    tp=sum(X_[cbind(li,lj)]>0)
    fn=sum(X)-tp#sum(!X_[cbind(li,lj)])
    fp=sum(X_@x>0)-tp
    cv=1-(fp+fn)/(tp+fn)
    stats=rbind(stats,cbind(tP,tp,fn,fp,cv,P=tp*1.0/(tp+fp),R=tp*1.0/(tp+fn)))
    }

   
    plot(stats[,'tP'],stats[,'R'],type='b',col='red',lwd=2,
    main=sprintf('topFiberM, dataset: \%s, 
         #Known facts:\%d','DBLP',sum(X)),ylab="",xlab='tP',
    xlim=c(0,1),ylim=c(0,1))
    HM=apply(stats,1,function(x){2/(1/x['P']+1/x['R'])})
    points(stats[,'tP'],stats[,'P'],col='blue',lwd=2,type='b')
    points(stats[,'tP'],HM,col='green',lwd=2,type='b')
    grid(nx=10, lty = "dotted", lwd = 2)
    legend(legend=c('Recall','Precision','Harmonic mean'),col=c('red','blue','green'),
    x=0.6,y=0.2,pch=1,cex=0.75,lwd=2)


}


}

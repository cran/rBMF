\name{topFiberM}
\alias{topFiberM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
topFiberM
}
\description{
 implements topFiberM Boolean matrix factorization algorithm. 
 topFiberM chooses in a greedy way the fibers (rows or columns) to represent the 
 entire matrix. Fibers are extended to rectangles according to a threshold on precision. 
 The search for these" top fibers" can continue beyond the required rank and according 
 to an optional parameter that defines the limit for this search.
}
\usage{
topFiberM(X, r = 2, tP = 0.5, verbose = 2, SR = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
the input boolean sparse matrix
}
  \item{r}{
rank (number of factors) required.
}
  \item{tP}{
  parameter to put threshold on precision
}
  \item{verbose}{
  integer value to control the appearance of messages.
  0 minimal messages will be showed. Default 2
}
  \item{SR}{
  search limit which defines the number iterations, minimum value is rank
  and maximum value is minimum number of columns and number of rows
}
}

\value{
List of the following four components:
\item{A}{Factor matrix A}
\item{B}{Factor matrix B}
\item{X1}{remaining uncovered ones, (False negatives)}
\item{tf}{dataframe logging of steps giving description of each factor,
contains index, based on column (2) / row (1), nnz, TP, FP}
}
\references{
Desouki, A. A., Roeder, M., & Ngomo, A. C. N. (2019). 
topFiberM: Scalable and Efficient Boolean Matrix Factorization.
arXiv preprint arXiv:1903.10326.
}
\author{
Abdelmoneim Amer Desouki
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{GreConD}}
 \code{\link{Asso_approximate}}
}
\examples{

data(DBLP)
    r=7
    tP=0.6
    X=DBLP
    Xb=X==1#Convert to boolean
    
 Res=topFiberM(Xb,r=r,tP=tP,SR=100,verbose=1)
    X_=Res$A \%*\% Res$B
    X_=as(X_,'TsparseMatrix')
    #Calculate metrics
    tempX=as(X,'TsparseMatrix')
    li=tempX@i[tempX@x==1]+1
    lj=tempX@j[tempX@x==1]+1
    tp=sum(X_[cbind(li,lj)]>0)
    fn=sum(X)-tp#sum(!X_[cbind(li,lj)])
    fp=sum(X_@x>0)-tp
    cv=1-(fp+fn)/(tp+fn)
    
print(sprintf("tp:\%d, fp:\%d,fn:\%d, Error:\%d, covered=\%.3f",tp,fp,fn,fn+fp,cv))
    
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ topFiberM }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{ Boolean Matrix Factorization }
% Use only one concept per line.

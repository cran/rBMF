\name{GreConDPlus}
\alias{GreConDPlus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GreConDPlus Boolean Matrix Factorization
}
\description{
implements GreConDPlus Boolean Matrix Factorization Algorithm.
}
\usage{
GreConDPlus(X, no_of_factors = NULL, w, verbose = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{X}{
input binary matrix to be factorized. 
}
  \item{no_of_factors}{
  number of factors of the result (optional).
}
  \item{w}{
  wieght parameter for the algorithm : represents the weight of type1 and 
  type2 errors
}
  \item{verbose}{
  integer value to control the appearance of messages.
  0 minimal messages will be showed. Default 2
}
}

\value{
List of the following four components:
\item{A}{Factor matrix A}
\item{B}{Factor matrix B}
}
\references{
Belohlavek R., Vychodil V.: 
 Discovery of optimal factors in binary data via a novel method of matrix decomposition. 
Journal of Computer and System Sciences 76(1)(2010), 3-20
}
\author{
Abdelmoneim Amer Desouki
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{topFiberM}}, 
 \code{\link{GreConD}},
 \code{\link{Asso_approximate}}
}

\examples{
\donttest{
data(DBLP)
X=DBLP
Xb=X==1
Res=GreConDPlus(Xb,7,w=4,verbose=1)

    X_=Res$A \%*\% Res$B
    X_=as(X_,'TsparseMatrix')
    X=as(X,'TsparseMatrix')
			li=X@i[X@x==1]+1
			lj=X@j[X@x==1]+1
			tp=sum(X_[cbind(li,lj)]>0)
			fn=sum(X)-tp#sum(!X_[cbind(li,lj)])
			fp=sum(X_@x>0)-tp
			cv=1-(fp+fn)/(tp+fn)
print(sprintf("tp:\%d, fp:\%d,fn:\%d, Error:\%d, covered=\%.3f",tp,fp,fn,fn+fp,cv))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
 \keyword{ GreConDPlus }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
 \concept{ Boolean Matrix Factorization }
% \concept{ ~cpt2 }
% Use only one concept per line.

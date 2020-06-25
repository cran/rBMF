\name{GreConD}
\alias{GreConD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GreConD Boolean matrix factorization
}
\description{
implements GreConD algorithm for Boolean matrix factorization.
returns A . B = X (if the no. of factors is not limited).
}
\usage{
GreConD(X, no_of_factors = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
input binary matrix to be factorized. 
}
  \item{no_of_factors}{
  number of factors of the result (optional).
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
 \code{\link{Asso_approximate}}
}

\examples{
data(DBLP)
X=DBLP
Xb=X==1
Res=GreConD(Xb,7)

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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
 \keyword{ GreConD }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{ Boolean Matrix Factorization }
% \concept{ ~cpt2 }
% Use only one concept per line.

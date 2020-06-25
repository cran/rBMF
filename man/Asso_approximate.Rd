\name{Asso_approximate}
\alias{Asso_approximate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Asso: Boolean Matrix Factorization
}
\description{
Given binary matrix S (m x n), and a scalar k (number of factors), Asso finds two matrices 
A (m x r), B (r x n) such taht S~= A * B 
}
\usage{
Asso_approximate(S, k, opti)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
input matrix to be factorized
}
  \item{k}{
number of factors (rank)
}
  \item{opti}{
options : list containing components verbose:control of showing messages , 
threshold:precision limit,penalty_overcovered,bonus_covered.
}
}

\value{
LIST of three components
  \item{B }{k x n factor matrix }
 \item{O }{n x k factor matrix}
 \item{D}{n x n matrix, the result from calculate_association}
%% ...
}
\references{
Miettinen, P., Mielikaeinen, T., Gionis, A., Das, G., & Mannila, H. (2008). The discrete basis problem. IEEE transactions on knowledge and data engineering, 20(10), 1348-1362.
}
\author{
Abdelmoneim Amer Desouki
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 See Also \code{\link{topFiberM}}
}
\examples{
data(DBLP)
X=DBLP
Xb=X==1
Res=Asso_approximate(Xb,7,list(threshold=0.5,penalty_overcovered=1,bonus_covered=1,verbose=0))

    X_=Res$O \%*\% Res$B
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
 \keyword{ Asso }
 \keyword{ DBP }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{ Boolean Matrix Factorization}
% \concept{ ~cpt2 }
% Use only one concept per line.

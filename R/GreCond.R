#  Copyright (C) 2020 Abdelmoneim Amer Desouki, 
#   Data Science Group, Paderborn University, Germany.
#  All right reserved.
#  Email: desouki@mail.upb.de
#
#  This file is part of rBMF package
#
#  rBMF is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  rBMF is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with rBMF.  If not, see <http://www.gnu.org/licenses/>.

#3/1/2019
## this code is converted from MATLAB code provided by the authors

 GreConD<-function( X, no_of_factors =NULL){#[ A, B ] 
# % GRECOND implements GreConD algorithm for Boolean matrix factorization 

# % usage: [A, B] = GreConD(X);
# % returns A \circ B = X (if the no. of factors is not limited)

# % if you are using this implementation please cite the following work
# % Belohlavek R., Vychodil V.: 
# % Discovery of optimal factors in binary data via a novel method of matrix decomposition. 
# % Journal of Computer and System Sciences 76(1)(2010), 3-20. 

 M = X==1#logical(X);
# [m, n] = size(M);
m=nrow(M)
n=ncol(M)
U = M;
k = 0;

# A = logical([]);
# B = logical([]);
A=Matrix::spMatrix(x=FALSE,i=1,j=1,nrow=m,ncol=1)
B=Matrix::spMatrix(x=FALSE,i=1,j=1,nrow=1,ncol=n)

while (sum(U)>0){#any(any(U))
    v = 0;
    d =rep(FALSE,n)#,false(1,n);
    d_old = d;
    d_mid = d#false(1,n);
    e = matrix(rep(TRUE,m),ncol=1)#true(m,1); #% extent for speed closure
    
    atr = which(colSums(U)>0); #% only not covered attributes
    
    while( 1==1){
         for (j in which(!d[atr])){
            if(!d[j]){#~d(j)
                # % computes the value of the cover function for the candidate factor
                # % inline function for speed
                # % arrow down (speed version)
                a = as.vector(e & M[,j]);
                # % arrow up
                # b = all(M(a,:),1);
                
                # if(sum(a)==0){
                    # b=rep(TRUE,n)#as MATLAB all, return all ones
                # }else{
                  b = colSums(M[a,,drop=FALSE])==sum(a);
                 # }
                # % coverage
                # cost = sum(sum(U(a,b)));
                cost = sum(U[a,b]);
                # %
                
                if (cost > v){
                    v = cost;
                    d_mid = b;
                    cc = a;
                }
            }
        }
        
        d = d_mid;
        e = cc;
        
        if (all(d==d_old)){
            break;
        }else{
            d_old = d;
        }
    }
    
    k = k + 1;
    print(k);
    # browser()
    # A = [A, cc];
    A[,k]=cc
    # B = [B; d];
    B[k,]=d;
    
    
    # % end if the no. of factors is reached
    if (!is.null(no_of_factors) && k==no_of_factors){
        break;
    }else{
        A=cbind(A,FALSE)
        B=rbind(B,FALSE)
    }
    
    # % delete already covered part
    U[cc, d] = FALSE;
}
return(list(A=A,B=B))
}
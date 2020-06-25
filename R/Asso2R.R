
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

## this code is converted from C code provided by the authors
# 19/6/2018

###
# con<-file('C:\\Users\\Abdelmonem\\Dropbox\\TF\\pauli\\mdl4bmf\\decompDB.txt','r')
read_matrix<-function(fname){
    con<-file(fname,'r')
    nr=as.integer(readLines(con,1))
    nc=as.integer(readLines(con,1))
    B=Matrix::sparseMatrix(x=FALSE,i=1,j=1,dims=c(nr,nc))
    for(r in 1:nr){
        ltxt=unlist(strsplit(readLines(con,1),split=' '))
        if(any(ltxt=='1'))
            B[r,ltxt=='1']=TRUE   
    }
    close(con)
return(B)
}

write_spmatrix<-function(A,fname){
    rs=rowSums(A)
    con<-file(fname,'w')
    print(fname)
    writeLines(as.character(nrow(A)),con)
    writeLines(as.character(ncol(A)),con)
    for(r in 1:length(rs)){
        writeLines(as.character(rs[r]),con)
        if(rs[r]>0)
            writeLines(paste(which(A[r,]==1),collapse=" "),con)
    }
    close(con)
}

#######
# sparse matrix is an array of linked lists (one per row)


calculate_association <- function(X,opti){
        cs=colSums(X)
        D=Matrix::sparseMatrix(x=FALSE,i=1,j=1,dims=c(ncol(X),ncol(X)))
       for(bit in 1:ncol(X)){ 
            if(opti$verbose>0) print(sprintf('col %d of %d ...',bit,ncol(X)))
            D[bit,]=(colSums(X[which(X[,bit]==1),,drop=FALSE])/cs[bit])>=opti$threshold
        }
        return(D)
}

solve_basis <-function(X,k,D,opti){
#D: is the result from calculate_association
#X~=O  o B
    nr=nrow(X)
    nc=ncol(X)
    B=Matrix::sparseMatrix(x=FALSE,i=1,j=1,dims=c(k,nc))
    O=Matrix::sparseMatrix(x=FALSE,i=1,j=1,dims=c(nr,k))#transposed
    
    covered=Matrix::sparseMatrix(x=FALSE,i=1,j=1,dims=c(nr,nc))#matrix(0,nr,nc)
    
    for (basis in 1: k) {
        if (opti$verbose > 0) {
          print(sprintf(" basis vector %d   ", basis));
        }
        best = -1;
        best_covers = 0;
        for (rw in 1:nc ) {# row in D:nc x nc
            element=which(D[rw,])
            rowcount=rowSums(X[,element,drop=FALSE]*opti$bonus_covered * (X[,element,drop=FALSE] - covered[,element,drop=FALSE]))- 
                     rowSums((!X[,element,drop=FALSE])*opti$penalty_overcovered * (! covered[,element,drop=FALSE]));
            covers =  sum(rowcount[rowcount>0])
          print(sprintf(" row:%d, covers=%d, bestcovers:%d",rw,covers,best_covers))
          if (covers > best_covers) {
            # /* we found the best */
               best_covers = covers;
            # /* 'best' is this row */
               best = rw;
               best_rowcount = rowcount
          }
        }#for rw
        
        if (best > -1) {
          # /* We have found the best - if best == -1, this basis vector will be
           # * empty.
           # */
           for(j in which(D[best,])){
            covered[cbind(which(best_rowcount>0),j)] = TRUE
           }            
          B[basis,]=D[best,]
        }
        O[cbind(which(best_rowcount > 0),basis)] = TRUE;
        # if (opti$verbose > 1) 
        print(paste0('basis:',basis,', best=',best,', best_covers=',best_covers))
    }
    return(list(B=B,O=O,D=D))
}
    
           # approximate
####------------------------------------
Asso_approximate<-function(S, k,opti){
  if (opti$verbose > 0) {
    print("Calculating associations...");
  }
  D = calculate_association(S, opti);

  if (is.null(D) )
    return (NULL);
    
  if (opti$verbose > 0) {
    print("Solving basis...");
  }
  return (solve_basis(S, k, D, opti));

}
######################
# calc_asso_rng <- function(X,opti){
#        
#         D=Matrix::sparseMatrix(x=FALSE,i=1,j=1,dims=c(ncol(X),ncol(X)))
#        for(bit in which(cs>0)){ 
#             if(opti$verbose>0) print(sprintf('col %d of %d ...',bit,ncol(X)))
#             ix=cs/cs[bit]<=2 & cs/cs[bit]>=opti$threshold
#             print(sum(ix))
#             D[bit,]=(colSums(X[which(X[,bit]==1),,drop=FALSE])/cs[bit])>=opti$threshold
#         }
#         return(D)
# }



  # if(basis==2){
                # rowcount1=rep(0,nr)            
                # Pos1=rep(0,nr)            
                # Neg1=rep(0,nr)       
                # Pos=rep(0,nr)            
                # Neg=rep(0,nr)       
                # covered1=read_matrix('C:\\Users\\Abdelmonem\\Dropbox\\TF\\pauli\\mdl4bmf\\covered1.txt')
                 
                # # covered=covered1==1
                # for (i in 1:nr) {
                    # #/* find best row */
                    # rowcount1[i] = sum((X[i,element])*opti$bonus_covered * (1 - covered[i,element])) - sum((!X[i,element])*opti$penalty_overcovered * (1 - covered[i,element]))  ;
                    # Pos[i]= opti$bonus_covered * sum(!covered[i,which(X[i,element])])
                    # Neg[i]= opti$penalty_overcovered * sum(!covered[i,which(!X[i,element])])
                    # Pos1[i]= opti$bonus_covered * sum(!covered1[i,which(X[i,element])])
                    # Neg1[i]= opti$penalty_overcovered * sum(!covered1[i,which(!X[i,element])])
                    # if (opti$verbose > 1) print(paste0('row=',rw,' i=',i,' rowcount=',rowcount[i]))
                    # #if (rowcount[i] > 0) covers = covers + rowcount[i];
                # }
                # browser()
            # }
          # covered1=read_matrix('C:\\Users\\Abdelmonem\\Dropbox\\TF\\pauli\\mdl4bmf\\covered1.txt')
          
          # con<-file('C:\\Users\\Abdelmonem\\Dropbox\\TF\\pauli\\mdl4bmf\\best_rowcount0.txt','r')
          # readLines(con,2)
          # bc0=as.integer(readLines(con))
          # close(con)
          # browser()
          # print(sprintf("SumCovered:%d",sum(covered)))
          
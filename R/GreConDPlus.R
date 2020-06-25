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


## this code is converted from MATLAB code provided by the authors
# 18/1/2019
 GreConDPlus<-function( X, no_of_factors=NULL, w ,verbose=2){#[ A, B ] = 

   cpp_uzaver <-function(D,y,U,M){
     
   };
   func <- 'List cpp_uzaver(Rcpp::NumericVector D, int y, NumericMatrix U, NumericMatrix M) {
int m=M.nrow();
int n=M.ncol();
int i,j;
int cost = 0;
//Rcout<<"n="<< n << ", m=" << m << "y=" << y << ", D[y-1]="<<D[y-1]<<std::endl;
D(y-1) = 1;
/*Rcout<<"U[0,0]"<< U[0,0] << ", U[0,1]" << U[0,1] << ",U[1,0]" << U[1,0] << std::endl;
for(i=0; i<m; i++){
  for(j=0; j<n; j++)  Rcout << U(i,j) << " ";
  Rcout<<std::endl;
  }
    for(j=0; j<n; j++) Rcout<<D(j)<<" ";
    */
Rcpp::NumericVector A(m);
Rcpp::NumericVector B(n);
// sipka dolu
  for(i=0; i<m; i++)
  {
    A(i) = 0;  
    for(j=0; j<n; j++) 
    {
    
     if (M(i,j) >= D(j)) A(i) = 1;
     else
     {
         A(i) = 0;
         break;
     }
    }
  } 

  for(j=0; j<n; j++)
  {
     B(j) = 0;      
     for(i=0; i<m; i++)
     {
        if (M(i,j) >= A(i)) B(j) = 1;
        else
        {
            B(j) = 0;
            break;
        }
     }
  }    
  
  // spocitame plochu 
  for(i=0; i<m; i++)
  {
    for(j=0; j<n; j++)
    {         
     if(A(i) && B(j) && U(i,j)) cost++;
    } 
  }  
  
  D(y-1) = 0;
  
  List ret;
ret["cost"] = cost;
ret["A"] = A;
ret["B"] = B;
return ret;
}'
   
   Rcpp::cppFunction(func)
   
M = X==1
# M = logical(vstup);
# [m, n] = size(M);
m=nrow(M)
n=ncol(M)
coverage = matrix(0,nrow=m, ncol=n);
U = M;
k = 0;

# a = logical([]);
# b = logical([]);
# E = logical([]); % expansion of the factor
# F = logical([]); % sloupcova expanze faktoru
a=NULL
b=NULL
E=NULL
FF=NULL
# % run the calculation
while(max(U)>0){
    max_cost = -Inf;
    D = rep(FALSE,n)#false(1,n);
    D_old = rep(FALSE,n)#false(1,n);
    flag = TRUE;#podminka
    D_between = rep(FALSE,n)#false(1,n);
    
    # % select the GreConD help factor
    while(flag){
        for(j in 1:n){
            if(!D[j]){
                # tmp = r_uzaver(D,j,U,M);
                tmp = cpp_uzaver(D,j,U+0,M+0);
  # if(k==1) browser();
                cost=tmp[[1]]; A=tmp[[2]]; B=tmp[[3]]
                if(verbose>2) print(sprintf("k=%d,j=%d,cost=%d, |A|=%d, U=%d",k,j,cost,sum(A),sum(U)));
                if (cost > max_cost){
                    max_cost = cost;
                    D_between = B;
                    C = A;
                }
            }
        }
        
        D = D_between;
        
        if(min(D==D_old)){
            flag = FALSE;
        }else{
            flag = TRUE;
            D_old = D;
        }
    }
    
    # % e, f represents the expansion factor factor [C, D]
    # browser()
    tmp = expansion(C, D, U, M, w);
    e=tmp[[1]]; f=tmp[[2]];
    C = C | e #or(C, e);
    D = D | f #or(D, f);
    
    if(is.null(a)){ 
        a = matrix(C,ncol=1)
    }else{
        a = cbind(a,C)#[a, C];
    }
    # b = [b; D];
    b=rbind(b,D)
    if(is.null(E)){ 
        E = matrix(e,ncol=1)
    }else{
        E = cbind(E,e)
    }#E = [E, e];
    # F = [F; f];
    FF = rbind(FF,f)
    
    k = k + 1;
    if(verbose>=2) print(k);
    
    # if k==pocet_faktoru && pocet_faktoru
        # break;
    # end
    if (!is.null(no_of_factors) && k==no_of_factors){#NB:no expansion for last factor!!!!!
        break;
    }
    # % remove from U and add to cover matrix
    ix=as.matrix(expand.grid(which(C==1),which(D==1)),ncol=2)
    U[ix] = FALSE;
    # browser()
    coverage[ix] = coverage[ix] + 1;
        
    # % overime whether it is possible to remove some factor (totally cancel)
    for (i in 1:nrow(b)){
        if (i >= nrow(b)){#size(b, 1)
            break;
        }
        
        ix2=as.matrix(expand.grid(which(a[,i]==1),which(b[i,]==1)),ncol=2)
        cc = coverage[ix2]#coverage(a(:,i), b(i,:));
        
        if (min(cc)>=2){#min(min(c(M(a(:,i), b(i,:))))) >= 2
            if(verbose>=2) print('taking the factor...');
            
            a=a[,-i]#a(:,i) = [];
            b=b[-i,]#b(i,:) = [];
            E=E[,-i]#E(:,i) = []; #% delete expansion
            FF=FF[-i,]#F(i,:) = [];
            # coverage(a(:,i), b(i,:)) = coverage(a(:,i), b(i,:)) - 1;
            coverage[ix2] = coverage[ix2] - 1
        }
    }
    
    # % verify whether the overlap error can be mitigated
    for(i in 1:nrow(b)){
        # % loop over all columns in expansion
        for( j in 1:n){
            if( FF[i,j]){
                z = rep(1,n)#false(1,n);
                z[j] = 1;
                ix3=as.matrix(expand.grid(which(a[,i]==1),which(z==1)),ncol=2)
                cc = coverage[ix3]#coverage(a(:,i), z);
                # % pokud je mozne odebrat expanzy (jednicky v expanzi jsou
                # % pokryty jinymi faktory)
                if(nrow(ix3)>0 && min(cc[M[ix3]]) >= 2){# min(min(c(M(a(:,i), z)))) >= 2
                    print('remove the expansion...');
                    b[i,j] = 0;
                    # coverage(a(:,i), z) = coverage(a(:,i), z) - 1;
                    coverage[ix3] = coverage[ix3] - 1;
                    FF[i,j] = 0;
                }
            }
        }
        
        # % loop over all lines in expansion
        for( j in 1:m){
            if (E[j,i]){
                z = rep(0,m)#false(m,1);
                z[j] = 1;
                # cc = coverage(z, b(i,:));
                ix4=as.matrix(expand.grid(which(z==1),which(b[i,]==1)),ncol=2)
                cc = coverage[ix4]#coverage(a(:,i), z);
                # browser()
                # % if it is possible to remove the expansions (they are one in expansion covered by other factors)
                if (nrow(ix4)>0 && min(cc[M[ix4]])>=2){#min(min(c(M(z, b(i,:))))) >= 2
                    print('removing the line from expansion...');
                    a[j,i] = 0;
                    # coverage(z, b(i,:)) = coverage(z, b(i,:)) - 1;
                    coverage[ix4] = coverage[ix4] - 1;
                    E[j,i] = 0;
                }
            }
        }
    }
      
}


A = a==1#logical(a);
B = b==1#logical(b);

return(list(A=A,B=B))
}



# % the function for expanding the concept by one row or column
 expansion <- function(A, B, U, M, w){#[E, F] 

# [m, n] = size(U);
m=nrow(U)
n=ncol(U)
E = rep(FALSE,m)#false(m,1);
FF = rep(FALSE,n)#false(1,n);

while( 1==1){    
    ix=as.matrix(expand.grid(which(A==1),which(B==1)),ncol=2)
    cover = sum(U[ix])#sum(sum(U(A, B)));
    overcover = sum(!M[ix])#sum(sum(~M(A, B)));
    cost = cover - w * overcover;
    co = 0; #% determines whether the column or row is a good one
    
    # % the price of all columns, except for those in B
    price = cover + colSums(U[A==1,]) - w * (overcover + colSums(!M[A==1,,drop=FALSE]))#sum(U(A,:)) - w * (overcover + sum(~M(A,:),1));
    price[B==1] = -Inf;
    colno = which.max(price);
   
    # % the price of the best column
    if (cost < price[colno]){
        co = 2;
    }
        
    # % the price of all lines except those in A
        # price = cover + sum(U(:,B),2) - w * (overcover + sum(~table(:,B),2));
    price = cover + rowSums(U[,B==1,drop=FALSE]) - w * (overcover + rowSums(!M[,B==1,drop=FALSE]));
    price[A==1] = -Inf;
    rowno = which.max(price);
   
    # % the price of the best line
    if (cost < price[rowno]){
        co = 1;
    }
    
   if (co == 0 ){#% zadna expanze
        break;
    }else{
        if (co == 1){# % expanze o radek
            A[rowno] = 1;
            E[rowno] = 1;
        }else{
            if( co == 2){# % expansion by column
                B[colno] = 1;
                FF[colno] = 1;
            }
        }
    }
}#while   
   return(list(E=E,FF=FF))
}

# library(Rcpp)
func <- 'List cpp_uzaver(Rcpp::NumericVector D, int y, NumericMatrix U, NumericMatrix M) {
int m=M.nrow();
int n=M.ncol();
int i,j;
int cost = 0;
//Rcout<<"n="<< n << ", m=" << m << "y=" << y << ", D[y-1]="<<D[y-1]<<std::endl;
D(y-1) = 1;
/*Rcout<<"U[0,0]"<< U[0,0] << ", U[0,1]" << U[0,1] << ",U[1,0]" << U[1,0] << std::endl;
for(i=0; i<m; i++){
  for(j=0; j<n; j++)  Rcout << U(i,j) << " ";
  Rcout<<std::endl;
  }
    for(j=0; j<n; j++) Rcout<<D(j)<<" ";
    */
Rcpp::NumericVector A(m);
Rcpp::NumericVector B(n);
// sipka dolu
  for(i=0; i<m; i++)
  {
    A(i) = 0;  
    for(j=0; j<n; j++) 
    {
    
     if (M(i,j) >= D(j)) A(i) = 1;
     else
     {
         A(i) = 0;
         break;
     }
    }
  } 

  for(j=0; j<n; j++)
  {
     B(j) = 0;      
     for(i=0; i<m; i++)
     {
        if (M(i,j) >= A(i)) B(j) = 1;
        else
        {
            B(j) = 0;
            break;
        }
     }
  }    
  
  // spocitame plochu 
  for(i=0; i<m; i++)
  {
    for(j=0; j<n; j++)
    {         
     if(A(i) && B(j) && U(i,j)) cost++;
    } 
  }  
  
  D(y-1) = 0;
  
  List ret;
ret["cost"] = cost;
ret["A"] = A;
ret["B"] = B;
return ret;
}'

# cppFunction(func)

# U=X==1
# D=rep(0,ncol(U))
# tmp=cpp_uzaver(D,1,U,U);

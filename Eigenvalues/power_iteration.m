%Given a diagonalizable matrix A, the algorithm will produce:
% -a number lambda, which is the greatest (in absolute value) eigenvalue of A, 
% -a nonzero vector x, which is a corresponding eigenvector of lambda
%The most time-consuming operation of the algorithm is the multiplication of matrix A by a vector, so it is effective for a very large sparse matrix, otherwise it may converge slowly.
function[k,r,x,lambda,flag]=power_iteration(A,x0,tol)
% k - number of iterations
% r - residual 
% x - eigen vector 
% lamda - eigen value
  k=0;
  x=x0;
  v=norm(A*x-(x'*A*x)*x);
  r(1)=norm(A*x-(x'*A*x)*x);
  flag=0;
  
  while  tol<v
      
      y=A*x;
      x=y/norm(y,2);
      v=norm(A*x-(x'*A*x)*x);
       k=k+1;
      r(k+1)=norm(A*x-(x'*A*x)*x);
      
      if k==150
          flag=1; 
        lambda=max(abs(eig(A)));
        break      
      end
  end

  c=A*x;
  if flag==0 
  lambda=c(1)/x(1);
  end

end 
%Successive over-relaxation (sor) is a method for solving linear systems A*x = b, where we use:
% tol - error tollerance 
% x0 - initial approximation (wil be rewritten and returned as a final approximation of x)
% k - number of steps needed for achieving good enough approximation (with respect to tol)
% w (omega) - relaxation parameter (used for lowering spectral radius of matrix A -> faster convergence from x0 to x)
% r - residual 

function [x0, k, r] = sor(A, b, tol, x0, w)

  k=0;
  v=norm(b,2);
  r(1)=norm(b-A*x0,2)/v;
  n=length(b);

  while r(k+1) > tol
      k=k+1;
      for i=1:n
          x0(i)=(1-w)*x0(i);
          pom=b(i);
          
          for j=1:i-1
              pom=pom-A(i,j)*x0(j);
          end
          
          for j=i+1:n
              pom=pom-A(i,j)*x0(j);
          end
          
          x0(i)=x0(i)+(pom*w)/A(i,i);
      end

  r(k+1)=norm(b-A*x0,2)/v;
end


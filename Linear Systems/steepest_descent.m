%Steepest descent algorithm can be used for solving A*x = b
%It's convergence can be very slow due to potential steps in directions which were already considered

function[x,k,v] = steepest_descent(A,b,tol,x0)

  r=b-A*x0;

  k=0;
  v(1)=norm(b-A*x0,2); %vector of residuals
  x=x0;

  while tol<(v/norm(b,2))
      
      k=k+1;
      g=r'*A*r;
      alpha=r'*r/g;
      
      x=x+alpha*r;
      r=r-alpha*A*r;
      v(k+1)=norm(b-A*x,2); 
      
  end

v=v/norm(b,2); 
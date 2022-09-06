function[x,k,res] = conjugate_gradient(A,b,x0,tol)
  
  k=0;
  d=b-A*x0;
  r=b-A*x0;
  res(1)=norm(r)/norm(b);
  
  while (norm(r)/norm(b))>tol
      
      temp=r.'*r;
      pom=A*d;
      alpha=temp/(d.'*pom);
      x0=x0+alpha*d;
      r=r-alpha*pom;
      beta=(r.'*r)/temp;
      d=r+beta*d;
      k=k+1;
      res(k+1)=norm(r)/norm(b);
  
  end

  x=x0;

end
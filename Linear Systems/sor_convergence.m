%This program is used to find the optimal omega (relaxation parameter for sor)
%rho refers to spectral radius of A, associated with given omega

function [rho_min,omega_opt] = sor_convergence(A)
  %From theory we know that sor converges for omega from [0,2]
  %I assumed a step increase of 0.01 but it can be altered to any desired value
  omega=0.01:0.01:2;
  n=length(omega);

  D=diag(diag(A));
  L=tril(A,-1);
  U=triu(A,1);
  rho_min=2000;
  
  for i=1:n;
      T=(D+omega(i)*L)\((1-omega(i))*D-omega(i)*U);
      rho(i)=max(abs(eig(T)));
          
      if (rho(i)<rho_min) 
          rho_min=rho(i);
          omega_opt=omega(i);    
      end
  end

  figure(1)
  plot(omega,rho)
  xlabel('omega')
  ylabel('rho')
  grid on

end


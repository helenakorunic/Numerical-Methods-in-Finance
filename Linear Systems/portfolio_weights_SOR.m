[n,~]=size(C);
e=ones(n,1);
x0=zeros(n,1);
tol=1e-5;
[~,omega_opt]=sor_convergence(C)
[x_sor,k_sor,r_sor]=sor(C,e,tol,x0,omega_opt)
weights=x_sor/(e'*x_sor)
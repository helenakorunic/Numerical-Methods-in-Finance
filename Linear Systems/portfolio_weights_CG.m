% C - covariance matrix
% mi - vector of returns for assets
% mi_P - expected portfolio return

[n,~]=size(C);
e=ones(n,1);
x0=zeros(n,1);
tol=1e-5;
mi_P=100;
[x_cg,k_cg,r_cg]=conjugate_gradient(C,e,tol,x0)
[y,~,~]=conjugate_gradient(C,mi,tol,x0);
a=e'*x;
b=e'*y;
c=mi'*y;
d=a*c-b^2;
g=(c*x-b*y)/d;
h=(a*y-b*x)/d;
weights=g+mi_P*h
%Inverse iteration algorithm is used to find an approximate eigenvector when an approximation to a corresponding eigenvalue is already known.
function [x,k,nor,lambda,conv]=inverse_iteration(A,x0,mu,tol)
    [n,~]=size(A);
    I=diag(diag(ones(n)));
    [L,U,P]=lu(A-mu*I);
    tmp=A*x0;
    k=0;
    nor(k+1)=norm(tmp-(x0'*tmp)*x0,2);
    while nor(k+1)>tol
        y=U\(L\(P*x0));
        x0=y/norm(y,2);
        k=k+1;
        t=A*x0;
        nor(k+1)=norm(t-(x0'*t)*x0,2);
    end
    plot(log(nor),'r*:');
    x=x0;
    %aproximation of an eigenvalue with Rayleigh coefficient
    lambda=(x'*A*x)/(x'*x);
    %speed of convergence
    ev=abs(eig(A)-mu*ones(n,1));
    ev=sort(esv,'ascend');
    conv=sv(1)/sv(2);
end
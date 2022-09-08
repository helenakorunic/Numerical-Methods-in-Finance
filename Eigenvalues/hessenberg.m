%Reduction to Hessenberg form algorithm
function [Q,H] = hessenberg(A)

n=size(A,1);
Q=eye(n);

if (n<=2)
    H=A;
    Q=eye(n);
else
    for i=1:n-2
        [v,beta,s]=gallery('house',A(i+1:n,i),0);
        A(i+1,i)=s;
        A(i+2:n,i)=zeros(n-i-1,1);
        for j=i+1:n 
            A(i+1:n,j)=A(i+1:n,j)-(beta*v'*A(i+1:n,j))*v;
        end
        for j=1:n 
            A(j,i+1:n)=A(j,i+1:n)-(beta*v'*A(j,i+1:n)')*v';
        end
        for j=2:n 
            Q(j,i+1:n)=Q(j,i+1:n)-(beta*v'*Q(j,i+1:n)')*v';
        end
        
    end
    H=A;
end
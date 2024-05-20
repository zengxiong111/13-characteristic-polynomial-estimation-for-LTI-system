function [A_hat B_hat] = A_B_least_square(U,X,N,m,T,p)

p=size(U,1);
n=size(X,1);
N=size(U,2);

U_hat = zeros(n+p,N-1);
X_hat = zeros(n+p,N-1);

for i =1:N-1
    U_hat(:,i) = [X(:,i);U(:,i)];
    X_hat(:,i) = [X(:,i+1);U(:,i+1)];
end

A_big_hat = X_hat*pinv(U_hat);
A_hat = A_big_hat(1:n,1:n);
B_hat = A_big_hat(1:n,n+1:n+p);


end
function p_hat = characteristic_polynomial_estimation_single(U,Y,n,num_traj)

m=size(Y,1);
p=size(U,1);

Y_hat  = zeros(m,num_traj);
U_hat  = zeros((m+p)*n,num_traj);

for i =1:num_traj
    U_hat_temp = [];
    %U_hat_temp  =  Y(:,(n+1)*(i-1)+1);
    for j=1:n
         U_hat_temp = [U_hat_temp;Y(:,(n+1)*(i-1)+j)];
    end
    for j=1:n
         U_hat_temp = [U_hat_temp;U(:,(n+1)*(i-1)+j)];
    end

    Y_hat(:,i) = Y(:,(n+1)*i);
    U_hat(:,i) = U_hat_temp;
end


p_big_hat = Y_hat*pinv(U_hat);
p_hat = p_big_hat(1:n);

  
end
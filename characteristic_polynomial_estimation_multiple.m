function p_hat = characteristic_polynomial_estimation_multiple(Y_m,U_m,num_traj,p,m,n)

 
Y_hat  = zeros(m,num_traj);
U_hat  = zeros((p+m)*n,num_traj);

for i =1:num_traj
    U_hat_temp = [];
     for j=1:n
         U_hat_temp = [U_hat_temp;Y_m(:,j,i)];
    end
    for j=1:n
         U_hat_temp = [U_hat_temp;U_m(:,j,i)];
    end
    U_hat(:,i) = U_hat_temp;
    
    Y_hat(:,i) = Y_m(:,n+1,i);
end


p_big_hat = Y_hat*pinv(U_hat);
p_hat = p_big_hat(1:n);

  
end
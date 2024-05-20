function q_hat = fabio_method(Y_m,U_m,num_traj,p,m,n)




Y_hat  = zeros(m*(n+1),num_traj);
U_hat  = zeros(p*n,num_traj);

for i =1:num_traj
    Y_hat_temp = [];
     for j=1: n+1
         Y_hat_temp = [Y_hat_temp;Y_m(:,j,i)];
     end
     Y_hat(:,i) = Y_hat_temp;
     
     U_hat_temp = [];
    for j=1:n
         U_hat_temp = [U_hat_temp;U_m(:,j,i)];
    end
    U_hat(:,i) = U_hat_temp;
end

G = zeros(1,n+1);
G(1,n+1) = 1;
H = [eye(n) zeros(n,1)];

N_d = 8*(n+n*m) + int16(16 * log(10));
U_d_hat  = zeros(p*n,N_d);

U_d_hat = U_hat(:,1:N_d);

X0_m = zeros(n,num_traj);
X0_m_d = zeros(n,N_d);

U_bar = [X0_m;U_hat];
U_bar_d = [X0_m_d;U_d_hat];

q_hat = -G * Y_hat * pinv(U_bar) * U_bar_d *...
    pinv([H*Y_hat*pinv(U_bar)*U_bar_d; U_d_hat ])*[eye(n);zeros(n*m,n)];

   
end
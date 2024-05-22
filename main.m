r=0.5;
n=3;
m=1; %outout dimension
p=1; %input dimension
T=n;
num_traj = 5000;
N = num_traj*(n+1);

X_1 = ones(n,1);

[A,B,C,D] = system_generation(r,m,n,p);

% noise level and input
sigma_u = 1;
sigma_w = 0.000;
sigma_z = 0.000;

Obs = C;
for i=1:n-1
    Obs = [Obs;C*A^i];
end

%The true characteristic polynomial parameters
%charpoly(A)
C*A^n*pinv(Obs)

%single trajectory estimation
[U_single,Y_single,X] = single_trajectory_generation(N,A,B,C,D,sigma_u,sigma_w,sigma_z,X_1);
p_hat_single = characteristic_polynomial_estimation_single(U_single,Y_single,n,num_traj) 

%Multiple trajectories estimation
Y_m = zeros(p,n+1,N);
U_m = zeros(m,n+1,N);
for i =1:N
    [U_temp,Y_temp,X] = single_trajectory_generation(n+2,A,B,C,D,sigma_u,sigma_w,sigma_z,X_1);
    Y_m(:,:,i) = Y_temp(:,1:n+1);
    U_m(:,:,i) = U_temp(:,1:n+1) ;
end

p_hat_multiple = characteristic_polynomial_estimation_multiple(Y_m,U_m,num_traj,p,m,n)

%Try the formular in Fabio's paper
q_hat = fabio_method(Y_m,U_m,num_traj,p,m,n)

% Estimation error
% fprintf('    The relative estimation error of G with ls:  %6.3E \n',norm(G_ls-G)./norm(G));
% fprintf('    The relative estimation error of A with ls:  %6.3E \n',norm(A_hat - A)./norm(A));


clear

r=0.5;
n=3;
m=1; %outout dimension
p=1; %input dimension
h = 2*n;
traj_length = 500;

% noise level and input
sigma_u = 1;
sigma_w = 0.001;
sigma_v = 0.001;

[A,B,C,D] = system_generation(r,m,n,p);

G = zeros(m,h*p);
G(:,1:p) = D;
for i = 1:h-1
    G(:,i*p + 1: (i+1)*p) = C*A^(i-1)*B;
end


%repeat the experiments to get the Monte Carlo estimation of the
%expectation
num_exp = 20;

x1_all = [1 1000];
num_x1 = size(x1_all,2);

traj_length_all = 1000*(1:10);
num_traj_length = size(traj_length_all,2);

all_loss_G_ls = zeros(num_exp,num_traj_length,num_x1);

mean_loss_G_ls = zeros(num_traj_length,num_x1);
var_loss_G_ls = zeros(num_traj_length,num_x1);


for i = 1:num_x1
    for j = 1:num_traj_length
        for k = 1:num_exp
            x_1 = x1_all(i);
            X_1 = x_1 * ones(n,1);
            
            %single trajectory estimation
            [U_single,Y_single,X] = single_trajectory_generation(traj_length_all(j)+h-1,A,B,C,D,...
                sigma_u,sigma_w,sigma_v,X_1);
             
            G_ls = G_least_square(U_single,Y_single,traj_length_all(j),m,h,p);
                   
            all_loss_G_ls(k,j,i) =  norm(G_ls - G, 'fro');
        end
        mean_loss_G_ls(j,i) = mean(all_loss_G_ls(:,j,i));
        var_loss_G_ls(j,i) = var(all_loss_G_ls(:,j,i));
    end
end

figure;
hold on;

errorbar(traj_length_all, mean_loss_G_ls(:,1),var_loss_G_ls(:,1),'-o','LineWidth',3); 
errorbar(traj_length_all,mean_loss_G_ls(:,2),var_loss_G_ls(:,2),'-o','LineWidth',3);  
 

legend('x1=1','x1=1000','Location','best' )
grid on;
ax = gca;
ax.LineWidth = 2;
ax.GridLineStyle = '--';
ax.GridAlpha = 0.8;
lgd.FontSize = 18;
xlabel('N','FontSize',18) ;
ylabel('Error of Markov Parameter Matrix','FontSize',18) ;
set(gca,'FontSize',20)
  

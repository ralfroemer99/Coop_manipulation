%%This script calculates the entire estimation model for the 2 agents conected
%%by some communication network, distributedly aggregating via the dynamic
%%average consensus. Before running, make sure LoadData contains the
%%correct input data and constrained_dynamics_init_2Manip_Model1 contains suitable
%%parameters. The final plots can show both, centralized consensus or
%%dynamic avg. consensus. Plots need possibly to be adapted.
clear
close all
clc
constrained_dynamics_init_2manip_Model1

% Number of timesteps
N = length(x_1);

%% Declare storage to save intermediate results
r_cons_save = zeros(N,6);
mo_cons_save = zeros(N,1);
jo_cons_save = zeros(N,6);

Sigma_r_cons_save = zeros(N,6);
Sigma_mo_cons_save = zeros(N,1);
Sigma_jo_cons_save = zeros(N,6);

r1_dyn_man1_save = zeros(N,3);
r1_dyn_man2_save = zeros(N,3);

Sigma_r1_dyn_man1_save = zeros(N,3);
Sigma_r1_dyn_man2_save = zeros(N,3);

r2_dyn_man1_save = zeros(N,3);
r2_dyn_man2_save = zeros(N,3);

Sigma_r2_dyn_man1_save = zeros(N,3);
Sigma_r2_dyn_man2_save = zeros(N,3);

mo_dyn_man1_save = zeros(N,1);
mo_dyn_man2_save = zeros(N,1);

Sigma_mo_dyn_man1_save = zeros(N,1);
Sigma_mo_dyn_man2_save = zeros(N,1);

jo_dyn_man1_save = zeros(N,6);
jo_dyn_man2_save = zeros(N,6);

Sigma_jo_dyn_man1_save = zeros(N,6);
Sigma_jo_dyn_man2_save = zeros(N,6);

% Storage for model accuracy calculation
t_man1_save = zeros(N, 3);
t_man2_save = zeros(N, 3);
PhiTheta_man1_save = zeros(N, 3);
PhiTheta_man2_save = zeros(N, 3);

to_man1_save = zeros(N, 3);
to_man2_save = zeros(N, 3);
PhiThetao_man1_save = zeros(N, 3);
PhiThetao_man2_save = zeros(N, 3);

% REMOVE: %%%%%%%%%
to2_man1_save = zeros(N, 3);
to2_man2_save = zeros(N, 3);

imp1_l_save = zeros(N, 6);
imp1_r_save = zeros(N, 6);
%%%%%%%%%%%%%%%%%%%

%% Starting values for the local estimators
% noise_init_gain = 0.5;
w_hat_man1 = [mu_0_1, mu_0_1, mu_0_1];
w_hat_man2 = [mu_0_2, mu_0_2, mu_0_2];

Sigma_k_man1 = repmat(Sigma_0_1,1,1,3);
Sigma_k_man2 = repmat(Sigma_0_1,1,1,3);

j1_hat_man1 = jo_hat_man1;
j2_hat_man1 = jo_hat_man1;
j3_hat_man1 = jo_hat_man1;

j1_hat_man2 = jo_hat_man2;
j2_hat_man2 = jo_hat_man2;
j3_hat_man2 = jo_hat_man2;

Sigmao_k_man1 = repmat(Sigma_j_0,1,1,3);
Sigmao_k_man2 = repmat(Sigma_j_0,1,1,3);

start_inertia = 300;

%% Iterate over dataset
for ii = 1:N
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1st ESTIMATOR (Translation) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Assign the posterior of the previous step as the prior of the current step
    Sigma_man1 = Sigma_k_man1;
    Sigma_man2 = Sigma_k_man2;
    
    % Calculate the input matrix Phi_s_mani and target vector target_mani. See eq.(13)
    [Phi_man1,t_man1] = model_trans(1,g,m,d,k,ddx_1(ii,:),dx_1(ii,:),x_1(ii,:),q1_init,ddx_1_des(ii,:),dx_1_des(ii,:),x_1_des(ii,:),ddx_2_des(ii,:),dx_2_des(ii,:),x_2_des(ii,:));
    [Phi_man2,t_man2] = model_trans(2,g,m,d,k,ddx_2(ii,:),dx_2(ii,:),x_2(ii,:),q2_init,ddx_2_des(ii,:),dx_2_des(ii,:),x_2_des(ii,:),ddx_1_des(ii,:),dx_1_des(ii,:),x_1_des(ii,:));
    
    % Store values for model accuracy calculation
    t_man1_save(ii, :) = t_man1;
    t_man2_save(ii, :) = t_man2;
    PhiTheta_man1_save(ii, :) = Phi_man1*theta_1_real;
    PhiTheta_man2_save(ii, :) = Phi_man2*theta_2_real;
    
    % Calculate the current paramter estimation vectors wa_hat, see eq. (18), (19)
    Sigma_k_man1(:,:,1) = Sigma_man1(:,:,1) - beta*(Sigma_man1(:,:,1)*(Phi_man1(1,:)'*Phi_man1(1,:))*Sigma_man1(:,:,1))/(1 + beta*Phi_man1(1,:)*Sigma_man1(:,:,1)*Phi_man1(1,:)');
    Sigma_k_man1(:,:,2) = Sigma_man1(:,:,2) - beta*(Sigma_man1(:,:,2)*(Phi_man1(2,:)'*Phi_man1(2,:))*Sigma_man1(:,:,2))/(1 + beta*Phi_man1(2,:)*Sigma_man1(:,:,2)*Phi_man1(2,:)');
    Sigma_k_man1(:,:,3) = Sigma_man1(:,:,3) - beta*(Sigma_man1(:,:,3)*(Phi_man1(3,:)'*Phi_man1(3,:))*Sigma_man1(:,:,3))/(1 + beta*Phi_man1(3,:)*Sigma_man1(:,:,3)*Phi_man1(3,:)');
    w_hat_man1(:,1) = Sigma_k_man1(:,:,1)*(Sigma_man1(:,:,1)^(-1)*w_hat_man1(:,1) + beta*Phi_man1(1,:)'*t_man1(1));
    w_hat_man1(:,2) = Sigma_k_man1(:,:,2)*(Sigma_man1(:,:,2)^(-1)*w_hat_man1(:,2) + beta*Phi_man1(2,:)'*t_man1(2));
    w_hat_man1(:,3) = Sigma_k_man1(:,:,3)*(Sigma_man1(:,:,3)^(-1)*w_hat_man1(:,3) + beta*Phi_man1(3,:)'*t_man1(3));
    
    Sigma_k_man2(:,:,1) = Sigma_man2(:,:,1) - beta*(Sigma_man2(:,:,1)*(Phi_man2(1,:)'*Phi_man2(1,:))*Sigma_man2(:,:,1))/(1 + beta*Phi_man2(1,:)*Sigma_man2(:,:,1)*Phi_man2(1,:)');
    Sigma_k_man2(:,:,2) = Sigma_man2(:,:,2) - beta*(Sigma_man2(:,:,2)*(Phi_man2(2,:)'*Phi_man2(2,:))*Sigma_man2(:,:,2))/(1 + beta*Phi_man2(2,:)*Sigma_man2(:,:,2)*Phi_man2(2,:)');
    Sigma_k_man2(:,:,3) = Sigma_man2(:,:,3) - beta*(Sigma_man2(:,:,3)*(Phi_man2(3,:)'*Phi_man2(3,:))*Sigma_man2(:,:,3))/(1 + beta*Phi_man2(3,:)*Sigma_man2(:,:,3)*Phi_man2(3,:)');
    w_hat_man2(:,1) = Sigma_k_man2(:,:,1)*(Sigma_man2(:,:,1)^(-1)*w_hat_man2(:,1) + beta*Phi_man2(1,:)'*t_man2(1));
    w_hat_man2(:,2) = Sigma_k_man2(:,:,2)*(Sigma_man2(:,:,2)^(-1)*w_hat_man2(:,2) + beta*Phi_man2(2,:)'*t_man2(2));
    w_hat_man2(:,3) = Sigma_k_man2(:,:,3)*(Sigma_man2(:,:,3)^(-1)*w_hat_man2(:,3) + beta*Phi_man2(3,:)'*t_man2(3));
    
    % Calculate the gPoE and ratio distributions, see eq. (23) to (30)
    if ii > 1
        w_cons_man1_old = w_cons_man1;
        Sigma_cons_man1_old = Sigma_cons_man1;
        w_cons_man2_old = w_cons_man2;
        Sigma_cons_man2_old = Sigma_cons_man2;    
    end
    [w_cons_man1, Sigma_cons_man1] = gPoE(w_hat_man1, Sigma_k_man1);
    [w_cons_man2, Sigma_cons_man2] = gPoE(w_hat_man2, Sigma_k_man2);
    if ii == 1
        w_cons_man1_old = w_cons_man1;
        Sigma_cons_man1_old = Sigma_cons_man1;
        w_cons_man2_old = w_cons_man2;
        Sigma_cons_man2_old = Sigma_cons_man2;    
    end
    
    % Save values of recent step for dyn avg cons
    r1_hat_old_man1 = r1_hat_man1;
    Sigma_r1_old_man1 = Sigma_r1_hat_man1;
    r2_hat_old_man1 = r2_hat_man1;
    Sigma_r2_old_man1 = Sigma_r2_hat_man1;
    mo_hat_old_man1 = mo_hat_man1;
    Sigma_mo_old_man1 = Sigma_mo_hat_man1;
    
    r1_hat_old_man2 = r1_hat_man2;
    Sigma_r1_old_man2 = Sigma_r1_hat_man2;
    r2_hat_old_man2 = r2_hat_man2;
    Sigma_r2_old_man2 = Sigma_r2_hat_man2;
    mo_hat_old_man2 = mo_hat_man2;
    Sigma_mo_old_man2 = Sigma_mo_hat_man2;
    
    psi_1_old = buildPsi([r1_hat_old_man1; r2_hat_old_man1; mo_hat_old_man1], [diag(Sigma_r1_old_man1); diag(Sigma_r2_old_man1); Sigma_mo_old_man1]);
    psi_2_old = buildPsi([r1_hat_old_man2; r2_hat_old_man2; mo_hat_old_man2], [diag(Sigma_r1_old_man2); diag(Sigma_r2_old_man2); Sigma_mo_old_man2]);

    % Calculate ratio distribution    
    [r1_hat_man1,r2_hat_man1,mo_hat_man1,Sigma_r1_hat_man1,Sigma_r2_hat_man1,Sigma_mo_hat_man1] = ratioDistribution(w_cons_man1, Sigma_cons_man1);
    [r2_hat_man2,r1_hat_man2,mo_hat_man2,Sigma_r2_hat_man2,Sigma_r1_hat_man2,Sigma_mo_hat_man2] = ratioDistribution(w_cons_man2, Sigma_cons_man2);
     
    psi_1 = buildPsi([r1_hat_man1; r2_hat_man1; mo_hat_man1], [diag(Sigma_r1_hat_man1); diag(Sigma_r2_hat_man1); Sigma_mo_hat_man1]);  
    psi_2 = buildPsi([r1_hat_man2; r2_hat_man2; mo_hat_man2], [diag(Sigma_r1_hat_man2); diag(Sigma_r2_hat_man2); Sigma_mo_hat_man2]);
    
    % Save previous consensus state
    xi_1_old = xi_1;
    xi_2_old = xi_2;
    
    % Dynamic average consensus, see eq. (32)
    xi_1 = dynAvgCons(1, xi_1_old, xi_2_old, psi_1, psi_1_old, AdjMat);
    xi_2 = dynAvgCons(2, xi_2_old, xi_1_old, psi_2, psi_2_old, AdjMat);
   
    % Obtain mean and variance of predictive distribution from consensus state, see eq. (34)
    [mu1_tilde, Sigma1_tilde] = ResolveChi(xi_1); 
    r1_dyn_man1_save(ii,:) = mu1_tilde(1:3);
    r2_dyn_man1_save(ii,:) = mu1_tilde(4:6);
    mo_dyn_man1_save(ii,:) = mu1_tilde(7);
    Sigma_r1_dyn_man1_save(ii,:) = Sigma1_tilde(1:3);
    Sigma_r2_dyn_man1_save(ii,:) = Sigma1_tilde(4:6);
    Sigma_mo_dyn_man1_save(ii,:) = Sigma1_tilde(7);
 
    [mu2_tilde, Sigma2_tilde] = ResolveChi(xi_2); 
    r1_dyn_man2_save(ii,:) = mu2_tilde(1:3);
    r2_dyn_man2_save(ii,:) = mu2_tilde(4:6);
    mo_dyn_man2_save(ii,:) = mu2_tilde(7);
    Sigma_r1_dyn_man2_save(ii,:) = Sigma2_tilde(1:3);
    Sigma_r2_dyn_man2_save(ii,:) = Sigma2_tilde(4:6);
    Sigma_mo_dyn_man2_save(ii,:) = Sigma2_tilde(7);
    
    % Centralized consensus
    [r1_hat, Sigma_r1_hat] = gPoE([r1_hat_man1,r1_hat_man2], cat(3, Sigma_r1_hat_man1,Sigma_r1_hat_man2));
    [r2_hat, Sigma_r2_hat] = gPoE([r2_hat_man1,r2_hat_man2], cat(3, Sigma_r2_hat_man1,Sigma_r2_hat_man2));
    [mo_hat, Sigma_mo_hat] = gPoE([mo_hat_man1,mo_hat_man2], [Sigma_mo_hat_man1,Sigma_mo_hat_man2]);
    
    r_cons_save(ii,:) = [r1_hat; r2_hat];
    mo_cons_save(ii,:) = mo_hat;
    
    Sigma_r_cons_save(ii,:) = [diag(Sigma_r1_hat)', diag(Sigma_r2_hat)'];
    Sigma_mo_cons_save(ii,:) = diag(Sigma_mo_hat);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2nd ESTIMATOR (Rotation) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     betao = gain(ii)*max(beta/(1+0.25*max(max([diag(Sigma_r1_hat),diag(Sigma_r2_hat)]))),.000001);
    % Get measurement variance
    if ii < start_inertia
        betao = 0;
    else
        betao = betao_init;
    end
    
    % Assign the posterior of the previous step as the prior of the current step
    Sigmao_man1 = Sigmao_k_man1;
    Sigmao_man2 = Sigmao_k_man2;
    
    % Calculate the input matrix Phi_s_mani and target vector target_man, see eq. (42)
%     [Phio_man1,to_man1] = model_rot(r1_hat,r2_hat,m,d,k,j_eef,delta,kappa,ddx_1(ii,:),dx_1(ii,:),x_1(ii,:),q1_init,q2_init,ddx_1_des(ii,:),dx_1_des(ii,:),x_1_des(ii,:),ddx_2_des(ii,:),dx_2_des(ii,:),x_2_des(ii,:));    
%     [Phio_man2,to_man2] = model_rot(r2_hat,r1_hat,m,d,k,j_eef,delta,kappa,ddx_2(ii,:),dx_2(ii,:),x_2(ii,:),q2_init,q1_init,ddx_2_des(ii,:),dx_2_des(ii,:),x_2_des(ii,:),ddx_1_des(ii,:),dx_1_des(ii,:),x_1_des(ii,:));
    [Phio_man1,to_man1,to2_man1,imp_l,imp_r] = model_rot(O_r_1,r2_hat,m,d,k,j_eef,delta,kappa,ddx_1(ii,:),dx_1(ii,:),x_1(ii,:),q1_init,q2_init,ddx_1_des(ii,:),dx_1_des(ii,:),x_1_des(ii,:),ddx_2_des(ii,:),dx_2_des(ii,:),x_2_des(ii,:), ddx_2(ii,:),dx_2(ii,:),x_2(ii,:),h_1(ii,:),h_2(ii,:));    
    [Phio_man2,to_man2,to2_man2] = model_rot(O_r_2,r1_hat,m,d,k,j_eef,delta,kappa,ddx_2(ii,:),dx_2(ii,:),x_2(ii,:),q2_init,q1_init,ddx_2_des(ii,:),dx_2_des(ii,:),x_2_des(ii,:),ddx_1_des(ii,:),dx_1_des(ii,:),x_1_des(ii,:), ddx_1(ii,:),dx_1(ii,:),x_1(ii,:),h_2(ii,:),h_1(ii,:));
    

    % REMOVE: Use real instead of estimated quantities
%     [Phio_man1,to_man1,to2_man1] = model_rot(O_r_1,O_r_2,m,d,k,j_eef,delta,kappa,ddx_1(ii,:),dx_1(ii,:),x_1(ii,:),q1_init,q2_init,ddx_1_des(ii,:),dx_1_des(ii,:),x_1_des(ii,:),ddx_2_des(ii,:),dx_2_des(ii,:),x_2_des(ii,:),ddx_2(ii,:),dx_2(ii,:),x_2(ii,:),h_1(ii,:),h_2(ii,:));    
%     [Phio_man2,to_man2,to2_man2] = model_rot(O_r_2,O_r_1,m,d,k,j_eef,delta,kappa,ddx_2(ii,:),dx_2(ii,:),x_2(ii,:),q2_init,q1_init,ddx_2_des(ii,:),dx_2_des(ii,:),x_2_des(ii,:),ddx_1_des(ii,:),dx_1_des(ii,:),x_1_des(ii,:),ddx_1(ii,:),dx_1(ii,:),x_1(ii,:),h_2(ii,:),h_1(ii,:));
    % REMOVE

    % Store values for model accuracy calculation
    to_man1_save(ii, :) = to_man1;
    to_man2_save(ii, :) = to_man2;
    PhiThetao_man1_save(ii, :) = Phio_man1*j_guess;
    PhiThetao_man2_save(ii, :) = Phio_man2*j_guess;
    
    % REMOVE: %%%%%%%
    to2_man1_save(ii, :) = to2_man1;
    to2_man2_save(ii, :) = to2_man2;
    imp1_l_save(ii,:) = imp_l;
    imp1_r_save(ii,:) = imp_r;
    %%%%%%%%%%%%%%%%%

    % Calculate the current paramter estimation vectors wa_hat, see eq. (18), (19)
    Sigmao_k_man1(:,:,1) = Sigmao_man1(:,:,1) - betao*(Sigmao_man1(:,:,1)*(Phio_man1(1,:)'*Phio_man1(1,:))*Sigmao_man1(:,:,1))/(1 + betao*Phio_man1(1,:)*Sigmao_man1(:,:,1)*Phio_man1(1,:)');
    Sigmao_k_man1(:,:,2) = Sigmao_man1(:,:,2) - betao*(Sigmao_man1(:,:,2)*(Phio_man1(2,:)'*Phio_man1(2,:))*Sigmao_man1(:,:,2))/(1 + betao*Phio_man1(2,:)*Sigmao_man1(:,:,2)*Phio_man1(2,:)');
    Sigmao_k_man1(:,:,3) = Sigmao_man1(:,:,3) - betao*(Sigmao_man1(:,:,3)*(Phio_man1(3,:)'*Phio_man1(3,:))*Sigmao_man1(:,:,3))/(1 + betao*Phio_man1(3,:)*Sigmao_man1(:,:,3)*Phio_man1(3,:)');
    j1_hat_man1 = Sigmao_k_man1(:,:,1)*(Sigmao_man1(:,:,1)^(-1)*j1_hat_man1 + betao*Phio_man1(1,:)'*to_man1(1));
    j2_hat_man1 = Sigmao_k_man1(:,:,2)*(Sigmao_man1(:,:,2)^(-1)*j2_hat_man1 + betao*Phio_man1(2,:)'*to_man1(2));
    j3_hat_man1 = Sigmao_k_man1(:,:,3)*(Sigmao_man1(:,:,3)^(-1)*j3_hat_man1 + betao*Phio_man1(3,:)'*to_man1(3));
    
    Sigmao_k_man2(:,:,1) = Sigmao_man2(:,:,1) - betao*(Sigmao_man2(:,:,1)*(Phio_man2(1,:)'*Phio_man2(1,:))*Sigmao_man2(:,:,1))/(1 + betao*Phio_man2(1,:)*Sigmao_man2(:,:,1)*Phio_man2(1,:)');
    Sigmao_k_man2(:,:,2) = Sigmao_man2(:,:,2) - betao*(Sigmao_man2(:,:,2)*(Phio_man2(2,:)'*Phio_man2(2,:))*Sigmao_man2(:,:,2))/(1 + betao*Phio_man2(2,:)*Sigmao_man2(:,:,2)*Phio_man2(2,:)');
    Sigmao_k_man2(:,:,3) = Sigmao_man2(:,:,3) - betao*(Sigmao_man2(:,:,3)*(Phio_man2(3,:)'*Phio_man2(3,:))*Sigmao_man2(:,:,3))/(1 + betao*Phio_man2(3,:)*Sigmao_man2(:,:,3)*Phio_man2(3,:)');
    j1_hat_man2 = Sigmao_k_man2(:,:,1)*(Sigmao_man2(:,:,1)^(-1)*j1_hat_man2 + betao*Phio_man2(1,:)'*to_man2(1));
    j2_hat_man2 = Sigmao_k_man2(:,:,2)*(Sigmao_man2(:,:,2)^(-1)*j2_hat_man2 + betao*Phio_man2(2,:)'*to_man2(2));
    j3_hat_man2 = Sigmao_k_man2(:,:,3)*(Sigmao_man2(:,:,3)^(-1)*j3_hat_man2 + betao*Phio_man2(3,:)'*to_man2(3));
    
    % Calculate the gPoE, see eq. (23), (24)
    jo_hat_old_man1 = jo_hat_man1;
    Sigma_jo_old_man1 = Sigma_jo_hat_man1;
    jo_hat_old_man2 = jo_hat_man2;
    Sigma_jo_old_man2 = Sigma_jo_hat_man2;
    
    [jo_hat_man1,Sigma_jo_hat_man1] = gPoE([j1_hat_man1,j2_hat_man1,j3_hat_man1], Sigmao_k_man1);    
    [jo_hat_man2,Sigma_jo_hat_man2] = gPoE([j1_hat_man2,j2_hat_man2,j3_hat_man2], Sigmao_k_man2); 
    
    % Dynamic average consensus
    xi_1_jo_old = xi_1_jo;
    xi_2_jo_old = xi_2_jo;
    
    xi_1_jo = dynAvgConsJo(1,xi_1_jo_old,xi_2_jo_old,jo_hat_man1,jo_hat_old_man1,Sigma_jo_hat_man1,Sigma_jo_old_man1,AdjMat);
    xi_2_jo = dynAvgConsJo(2,xi_2_jo_old,xi_1_jo_old,jo_hat_man2,jo_hat_old_man2,Sigma_jo_hat_man2,Sigma_jo_old_man2,AdjMat);
    
    [jo_dyn_man1_save(ii,:),Sigma_jo_dyn_man1_save(ii,:)] = ResolveChiJo(xi_1_jo);
    [jo_dyn_man2_save(ii,:),Sigma_jo_dyn_man2_save(ii,:)] = ResolveChiJo(xi_2_jo);
    
    % Centralized consensus
    [jo_hat, Sigma_jo_hat] = gPoE([jo_hat_man1,jo_hat_man2], cat(3,Sigma_jo_hat_man1,Sigma_jo_hat_man2));
    jo_cons_save(ii,:) = jo_hat;
    Sigma_jo_cons_save(ii,:) = diag(Sigma_jo_hat);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Print progress
    if mod(ii,1000) == 0
        disp(ii + " steps done")
    end
end

% Plot results
plot_results

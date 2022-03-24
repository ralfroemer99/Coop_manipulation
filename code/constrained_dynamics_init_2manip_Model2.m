%%This script initializes the entire estimation of the decoupled estimator.
%%Before starting, make sure you want to use the preset noise parameter, otherwise
%%noise has to be adapted. Impedance controller can be adapted to suit the
%%data's origin. Initial states of the dyn. avg. cons. are adaptive.
%%Prepares for the execution of the script 'OfflineEstimatorModel2'.
m_o = 2.493;                %object mass in kg
g = [0;0;-9.81];            %gravity vector

beta = .05;          %variance assumed for estimator

% grasping points
O_r_1 = [0.01;   -0.24;  0.065]; %c_1 location w.r.t {o}
O_r_2 = [0.01;  0.32;  0.06]; %c_2 location w.r.t {o}
j_guess = [0.1;0.001;0.001;0.05;0.001;0.1];

Sel = eye(31);
Sel([8,11,12,26,29,30],:) = [];

%Starting values for the different manipulators' first estimator
mu_0_1 = Sel*[m_o;m_o*O_r_1;m_o*kron(O_r_1,O_r_1);j_guess;O_r_2-O_r_1;kron(O_r_2-O_r_1,O_r_2-O_r_1)];  
mu_0_2 = Sel*[m_o;m_o*O_r_2;m_o*kron(O_r_2,O_r_2);j_guess;O_r_1-O_r_2;kron(O_r_1-O_r_2,O_r_1-O_r_2)];
mu_0_1_opt = Sel*[m_o;m_o*O_r_1;m_o*kron(O_r_1,O_r_1);j_guess;O_r_2-O_r_1;kron(O_r_2-O_r_1,O_r_2-O_r_1)];

mu_0_2_opt = Sel*[m_o;m_o*O_r_2;m_o*kron(O_r_2,O_r_2);j_guess;O_r_1-O_r_2;kron(O_r_1-O_r_2,O_r_1-O_r_2)];
%Starting covariance for first estimator
Sigma_0_1 = blkdiag(1.5,diag(0.35*ones(9,1)),diag(0.05*ones(6,1))+0.0001,diag(0.3*ones(9,1)));

%Impedance controller
m = 10;             %mass of eefs in kg
j_eef = 0.5;        %inertia of eefs
d = 150;            %translational damping of eefs
delta = 1;        %rotational damping of eefs
k = 500;            %translational stiffness of eefs
kappa = 0.14;       %rotational stiffness of eefs

M = [m*eye(3),  zeros(3);
     zeros(3),  j_eef*eye(3)];

D = [d*eye(3),  zeros(3);
     zeros(3),  delta*eye(3)]; 
 
 %Cut data correctly
 ExtractFromStates
 ddx_1 = ddx_left(:,flag_coop:flag_coop_done)';
 dx_1 = dx_left(:,flag_coop:flag_coop_done)';
 x_1 = x_left(:,flag_coop:flag_coop_done)';
 
 ddx_2 = ddx_right(:,flag_coop:flag_coop_done)';
 dx_2 = dx_right(:,flag_coop:flag_coop_done)';
 x_2 = x_right(:,flag_coop:flag_coop_done)';
 
 ddx_jd = ddx_jd(:,flag_coop:flag_coop_done)';
 dx_jd = dx_jd(:,flag_coop:flag_coop_done)';
 x_jd = x_jd(:,flag_coop:flag_coop_done)';
 
 
 ddx_1_desired = ddx_jd(:,1:6);
 dx_1_desired = dx_jd(:,1:6);
 x_1_desired = x_jd(:,1:7);
 
 ddx_2_desired = ddx_jd(:,7:12);
 dx_2_desired = dx_jd(:,7:12);
 x_2_desired = x_jd(:,8:14);
 
 %Initial quaternions
 q1_initial = x_1(1,4:7);
 q2_initial = x_2(1,4:7);
 
 %AdjMat
 AdjMat = [1/2, 1/2;
           1/2, 1/2];
 
 
 % %Starting values for dynamic average consensus EEF. 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m1 = mu_0_1(2)/mu_0_1(1);
s1 = m1^2*((Sigma_0_1(2,2)/mu_0_1(2)^2)+(Sigma_0_1(1,1)/mu_0_1(1)^2));
a1 = m1/s1;
m2 = mu_0_1(3)/mu_0_1(1);
s2 = m2^2*((Sigma_0_1(3,3)/mu_0_1(3)^2)+(Sigma_0_1(1,1)/mu_0_1(1)^2));
a2 = m2/s2;
m3 = mu_0_1(4)/mu_0_1(1);
s3 = m3^2*((Sigma_0_1(4,4)/mu_0_1(4)^2)+(Sigma_0_1(1,1)/mu_0_1(1)^2));
a3 = m3/s3;

r1_hat_man1 = [m1;m2;m3];
xi_1_r1 = [a1;a2;a3;1/s1;1/s2;1/s3];
Sigma_r1_hat_man1 = diag([s1,s2,s3]);

m1 = mu_0_1(17) + mu_0_1(2)/mu_0_1(1);
s1 = Sigma_0_1(17,17) + m1^2*((Sigma_0_1(2,2)/mu_0_1(2)^2)+(Sigma_0_1(1,1)/mu_0_1(1)^2));
a1 = m1/s1;
m2 = mu_0_1(18) + mu_0_1(3)/mu_0_1(1);
s2 = Sigma_0_1(18,18) + m2^2*((Sigma_0_1(3,3)/mu_0_1(3)^2)+(Sigma_0_1(1,1)/mu_0_1(1)^2));
a2 = m2/s2;
m3 = mu_0_1(19) + mu_0_1(4)/mu_0_1(1);
s3 = Sigma_0_1(19,19) + m3^2*((Sigma_0_1(4,4)/mu_0_1(4)^2)+(Sigma_0_1(1,1)/mu_0_1(1)^2));
a3 = m3/s3;

r2_hat_man1 = [m1;m2;m3];
xi_1_r2 = [a1;a2;a3;1/s1;1/s2;1/s3];
Sigma_r2_hat_man1 = diag([s1,s2,s3]);

jo_hat_man1 = mu_0_1(11:16);
Sigma_jo_hat_man1 = Sigma_0_1(11:16,11:16)+0.0001;
xi_1_jo = [Sigma_jo_hat_man1^(-1)*mu_0_1(11:16);1./(reshape(Sigma_jo_hat_man1,[36,1]))];


mo_hat_man1 = mu_0_1(1);
xi_1_mo = [mu_0_1(1)/(diag(Sigma_0_1(1,1)));1./(Sigma_0_1(1,1))];
Sigma_mo_hat_man1 = Sigma_0_1(1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Starting values for dynamic average consensus EEF. 2
m1 = mu_0_2(17) + mu_0_2(2)/mu_0_2(1);
s1 = Sigma_0_1(17,17) + m1^2*((Sigma_0_1(2,2)/mu_0_2(2)^2)+(Sigma_0_1(1,1)/mu_0_2(1)^2));
a1 = m1/s1;
m2 = mu_0_2(18) + mu_0_2(3)/mu_0_2(1);
s2 = Sigma_0_1(18,18) + m2^2*((Sigma_0_1(3,3)/mu_0_2(3)^2)+(Sigma_0_1(1,1)/mu_0_2(1)^2));
a2 = m2/s2;
m3 = mu_0_2(19) + mu_0_2(4)/mu_0_2(1);
s3 = Sigma_0_1(19,19) + m3^2*((Sigma_0_1(4,4)/mu_0_2(4)^2)+(Sigma_0_1(1,1)/mu_0_2(1)^2));
a3 = m3/s3;

r1_hat_man2 = [m1;m2;m3];
xi_2_r1 = [a1;a2;a3;1/s1;1/s2;1/s3];
Sigma_r1_hat_man2 = diag([s1,s2,s3]);

m1 = mu_0_2(2)/mu_0_2(1);
s1 = m1^2*((Sigma_0_1(2,2)/mu_0_2(2)^2)+(Sigma_0_1(1,1)/mu_0_2(1)^2));
a1 = m1/s1;
m2 = mu_0_2(3)/mu_0_2(1);
s2 = m2^2*((Sigma_0_1(3,3)/mu_0_2(3)^2)+(Sigma_0_1(1,1)/mu_0_2(1)^2));
a2 = m2/s2;
m3 = mu_0_2(4)/mu_0_2(1);
s3 = m3^2*((Sigma_0_1(4,4)/mu_0_2(4)^2)+(Sigma_0_1(1,1)/mu_0_2(1)^2));
a3 = m3/s3;

r2_hat_man2 = [m1;m2;m3];
xi_2_r2 = [a1;a2;a3;1/s1;1/s2;1/s3];
Sigma_r2_hat_man2 = diag([s1,s2,s3]);

jo_hat_man2 = mu_0_2(11:16);
Sigma_jo_hat_man2 = Sigma_0_1(11:16,11:16)+0.0001;
xi_2_jo = [Sigma_jo_hat_man2^(-1)*mu_0_2(11:16);1./(reshape(Sigma_jo_hat_man2,[36,1]))];

mo_hat_man2 = mu_0_2(1);
xi_2_mo = [mu_0_2(1)/(diag(Sigma_0_1(1,1)));1./(Sigma_0_1(1,1))];
Sigma_mo_hat_man2 = Sigma_0_1(1,1);

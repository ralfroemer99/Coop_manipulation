% This script initializes the entire estimation of the decoupled estimator.
% Before starting, make sure you want to use the preset noise parameter, otherwise
% noise has to be adapted. Impedance controller can be adapted to suit the
% data's origin. Initial states of the dyn. avg. cons. are adaptive.
% Prepares for the execution of the script 'OfflineEstimatorModel1'.
%% Define parameters
% Sampling time
TA = 0.01;

% Object mass
m_o = 3.6+2*0.4;        % object plus flanges and screws
m_o = m_o + 6*1.043;    % Additional weights

% Gravity vector
g = [0;0;-9.81];

% Variance assumed for the estimator. The larger the betas, the more the measurements are weighted.
beta = 0.05;
betao_init = 0.05;
beta_sim = 0;

% Impedance controller parameters
m = 1;                              % mass of eefs in kg
k = 1000;                           % translational stiffness of eefs
d = 0.7*2*sqrt(m*k);                % translational damping of eefs
j_eef = 1;                          % inertia of eefs
kappa = 100;                        % rotational stiffness of eefs
delta = 0.7*2*sqrt(j_eef*kappa);    % rotational damping of eefs

% Consensus before ratio distribution
cons_before_rd = 0;

% Plot simulink data
plot_measured_data = 0;
plot_tracking_errors = 1;
plot_man_idx = 0;

% grasping points
O_r_1 = [0; -0.3225;  0]; %c_1 location w.r.t {o}
O_r_2 = [0;  0.3225;  0]; %c_2 location w.r.t {o}
% O_r_1 = [0; -0.35;  0]; %c_1 location w.r.t {o}
% O_r_2 = [0;  0.35;  0]; %c_2 location w.r.t {o}

% Initial estimate of the object inertia tensor
j_guess = guessInertia(4.4);

%% Starting values for the first estimator
mo_init = 2; % 2
O_r_1_init = O_r_1+[-0.11; 0.29; 0.16]; % [-0.11; 0.29; 0.16]
O_r_2_init = O_r_2+[-0.18; 0.25; 0.03]; % [-0.18; 0.25; 0.03]

mu_0_1 = [O_r_2_init-O_r_1_init;mo_init*O_r_1_init;mo_init];
mu_0_2 = [O_r_1_init-O_r_2_init;mo_init*O_r_2_init;mo_init];

theta_1_real = [O_r_2-O_r_1;m_o*O_r_1;m_o];
theta_2_real = [O_r_1-O_r_2;m_o*O_r_2;m_o];

% Starting covariance for first estimator
Sigma_0_1 = blkdiag(diag(0.35*ones(3,1)),diag(0.3*ones(3,1)),1.5);
% Sigma_0_1 = diag([0.35*ones(3,1); 0.3*ones(3,1); 1.5]);
Sigma_j_0 = diag(0.05*ones(6,1)) + 0.00000001;

%% Cut data correctly
% Extract from Simulink data
extract_from_simulink
ddx_1 = ddx_left(:,flag_coop:flag_coop_done)';
ddx_1_filt = ddx_left_filt(:,flag_coop:flag_coop_done)';
dx_1 = dx_left(:,flag_coop:flag_coop_done)';
dx_1_filt = dx_left_filt(:,flag_coop:flag_coop_done)';
x_1 = x_left(:,flag_coop:flag_coop_done)';
x_1_filt = x_left_filt(:,flag_coop:flag_coop_done)';

ddx_2 = ddx_right(:,flag_coop:flag_coop_done)';
ddx_2_filt = ddx_right_filt(:,flag_coop:flag_coop_done)';
dx_2 = dx_right(:,flag_coop:flag_coop_done)';
dx_2_filt = dx_right_filt(:,flag_coop:flag_coop_done)';
x_2 = x_right(:,flag_coop:flag_coop_done)';
x_2_filt = x_right_filt(:,flag_coop:flag_coop_done)';

h_1 = h_left(:,flag_coop:flag_coop_done)';
h_2 = h_right(:,flag_coop:flag_coop_done)';

% Plot measured and filtered data
t = 0:0.01:(0.01*(size(x_1,1)-1));
if plot_measured_data
    plot_data(x_1, dx_1, ddx_1, x_1_filt, dx_1_filt, ddx_1_filt, [1, 2, 3], t, 1);
    plot_data(x_2, dx_2, ddx_2, x_2_filt, dx_2_filt, ddx_2_filt, [1, 2, 3], t, 2);
end

% Use filtered data in the sequel
x_1 = x_1_filt(:,1:7);
dx_1 = dx_1_filt(:,1:6);
ddx_1 = ddx_1_filt(:,1:6);
x_2 = x_2_filt(:,1:7);
dx_2 = dx_2_filt(:,1:6);
ddx_2 = ddx_2_filt(:,1:6);

ddx_jd = ddx_jd(:,flag_coop:flag_coop_done)';
dx_jd = dx_jd(:,flag_coop:flag_coop_done)';
x_jd = x_jd(:,flag_coop:flag_coop_done)';

ddx_1_des = ddx_jd(:,1:6);
dx_1_des = dx_jd(:,1:6);
x_1_des = x_jd(:,1:7);

ddx_2_des = ddx_jd(:,7:12);
dx_2_des = dx_jd(:,7:12);
x_2_des = x_jd(:,8:14);

if plot_tracking_errors
    plot_tracking_error(x_1, x_1_des, t);
    plot_tracking_error(x_2, x_2_des, t);
end

if plot_man_idx
    tmp = size(man_idx1, 2);
    figure();
    plot(1:1:tmp, man_idx1); hold on
    plot(1:1:tmp, man_idx2);
end

%Initial quaternions
q1_init = x_1(1,4:7);
q2_init = x_2(1,4:7);

%AdjMat
AdjMat = [1/2, 1/2;
    1/2, 1/2];

%% Compute starting values for consensus states
%%%% EEF 1 %%%%
[r1_hat_man1, r2_hat_man1, mo_hat_man1, Sigma_r1_hat_man1, Sigma_r2_hat_man1, Sigma_mo_hat_man1] = ratioDistribution(mu_0_1, Sigma_0_1);
xi_1 = buildPsi([r1_hat_man1; r2_hat_man1; mo_hat_man1], [diag(Sigma_r1_hat_man1); diag(Sigma_r2_hat_man1); diag(Sigma_mo_hat_man1)]);
    
%%%% EEF 2 %%%%  
[r2_hat_man2, r1_hat_man2, mo_hat_man2, Sigma_r1_hat_man2, Sigma_r2_hat_man2, Sigma_mo_hat_man2] = ratioDistribution(mu_0_2, Sigma_0_1);
xi_2 = buildPsi([r1_hat_man2; r2_hat_man2; mo_hat_man2], [diag(Sigma_r1_hat_man2); diag(Sigma_r2_hat_man2); diag(Sigma_mo_hat_man2)]);

jo_hat_man1 = j_guess + [-0.1; -0.03; 0.02; 0.08; 0.03; 0.12];
Sigma_jo_hat_man1 = Sigma_j_0;
xi_1_jo = [Sigma_jo_hat_man1^(-1)*jo_hat_man1;1./(reshape(Sigma_jo_hat_man1,[36,1]))];

jo_hat_man2 = j_guess + [-0.13; -0.02; 0.05; 0.09; 0.04; 0.14];
Sigma_jo_hat_man2 = Sigma_j_0;
xi_2_jo = [Sigma_jo_hat_man1^(-1)*jo_hat_man2;1./(reshape(Sigma_jo_hat_man2,[36,1]))];

%% Plot data from simulink
function plot_data(x_1, dx_1, ddx_1, x_1_filt, dx_1_filt, ddx_1_filt, states, t, man)
n_states = length(states);

% Plot position
for ii = 1:n_states
    figure(3*n_states*(man-1) + ii)
    plot(t,x_1(:,states(ii)),'Color','red'); hold on
    plot(t,x_1_filt(:,states(ii)),'Color','blue');
    switch states(ii)
        case 1
            title("x-position (manipulator " + man + ")")
        case 2
            title("y-position (manipulator " + man + ")")
        case 3
            title("z-position (manipulator " + man + ")")
    end
end

% Plot velocities
n_filt = size(dx_1_filt,2)/6;
for ii = 1:n_states
    figure(3*n_states*(man-1) + ii + n_states)
    plot(t,dx_1(:,states(ii)),'Color','red'); hold on
    for jj = 1:n_filt
        plot(t,dx_1_filt(:,states(ii)+(jj-1)*6),'Color','blue'); hold on
    end
    switch states(ii)
        case 1
            title("x-velocity (manipulator " + man + ")")
        case 2
            title("y-velocity (manipulator " + man + ")")
        case 3
            title("z-velocity (manipulator " + man + ")")
    end
end

% Plot accelerations
figure(3+3*(man-1))
for ii = 1:n_states
    figure(3*n_states*(man-1) + ii + 2*n_states)
    plot(t,ddx_1(:,states(ii)),'Color','red'); hold on
    for jj = 1:n_filt
        plot(t,ddx_1_filt(:,states(ii)+(jj-1)*6),'Color','blue'); hold on
    end
    switch states(ii)
        case 1
            title("x-acceleration (manipulator " + man + ")")
        case 2
            title("y-acceleration (manipulator " + man + ")")
        case 3
            title("z-acceleration (manipulator " + man + ")")
    end
end

end

function plot_tracking_error(x, x_des, t)
figure()
plot(t, x(:,1), 'Color', 'red'); hold on
plot(t, x_des(:,1), 'Color', 'blue');
legend('x', 'x_d')
title('x-position')
figure()
plot(t, x(:,2), 'Color', 'red'); hold on
plot(t, x_des(:,2), 'Color', 'blue');
legend('y', 'y_d')
title('y-position')
figure()
plot(t, x(:,3), 'Color', 'red'); hold on
plot(t, x_des(:,3), 'Color', 'blue');
legend('z', 'z_d')
title('z-position')
figure()
tmp = quat2axang(x(:,4:7));
angle = rad2deg(tmp(:,4));
tmp = quat2axang(x_des(:,4:7));
angle_des = rad2deg(tmp(:,4));
plot(angle,'g'); hold on
plot(angle_des,'g--');
legend('Angle','Desired angle');
end

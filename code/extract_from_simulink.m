% This script extracts the necessary data from the variables saved during
% manipulation. 
% Load data
x = load('./data/test20/x.mat');
x_d = load('./data/test20/x_des.mat');
h = load('./data/test20/h.mat');
% q = load('./data/test15/q.mat');

% Downsample to 100 Hz
x.x_imp = x.x_imp(:,1:10:end);
x_d.h = x_d.h(:,1:10:end);
h.h = h.h(:,1:10:end);

% man_idx1 = get_manipulation_index(q.x_imp(2:8,:));
% man_idx2 = get_manipulation_index(q.x_imp(9:15,:));

%% Extract states
% Left robot
x_left = x.x_imp(2:8,:);
dx_left = get_dx(x_left,TA);
ddx_left = get_ddx(dx_left,TA);

% Apply moving mean filtering
% k = [31, 51, 101];
window_size = 7;
x_left_filt = movmean(x_left, 3, 2);
x_left_filt(4:7,:) = x_left_filt(4:7,:)./VecNorm(x_left_filt(4:7,:),2,1);
[dx_left_filt, ddx_left_filt] = derive_and_mean_filt(x_left_filt, window_size, TA);

% Right robot
x_right = x.x_imp(9:15,:);
dx_right = get_dx(x_right,TA);
ddx_right = get_ddx(dx_right,TA);

% Apply moving mean filtering
x_right_filt = movmean(x_right, 3, 2);
[dx_right_filt, ddx_right_filt] = derive_and_mean_filt(x_right_filt, window_size, TA);

% Desired values
x_jd = x_d.h(2:15,:);
x_jd(1:3,:) = 0.001*x_jd(1:3,:);        % Convert distances from mm into m.
x_jd(8:10,:) = 0.001*x_jd(8:10,:);
dx_jd = [get_dx(x_jd(1:7,:),TA); get_dx(x_jd(8:14,:),TA)];
ddx_jd = [get_ddx(dx_jd(1:6,:),TA); get_ddx(dx_jd(7:12,:),TA)];

% Wrenches
h_left = h.h(2:7,:);
h_right = h.h(8:13,:);
h_left = movmean(h_left,5,2);
h_right = movmean(h_right,5,2);

% Detect beginning and end of manipulation
flag_coop = find(dx_jd(4,100:end)~=0,1)-100;
% flag_coop = 10000;
% flag_coop = find(dx_jd(4,100:end)~=0,1);
disp('Starting of dynamic cooperative manipulation at')
disp(['Sample: ',num2str(flag_coop)])
% TODO: find period of zero desired motion longer than waiting time between the different motion phases.
flag_coop_done = size(x_left,2); 
disp('Ending dynamic cooperative manipulation at')
disp(['Sample: ',num2str(flag_coop_done)])

function dx = get_dx(x,TA)
% Compute the linear and angular velocity from the pose. The rotation is
% expressed in quaternions.
    N = size(x,2);
    dx = zeros(6,N);
    % Get linear velocity by computing v(t) = (x(t+1)-x(t-1))/(2*TA)
    for ii = 2:N-1
        for kk = 1:3
            dx(kk,ii) = (x(kk,ii+1)-x(kk,ii-1))/(2*TA);
        end
    end
    
    % Get angular velocity by computing w(t) =
    % (2*(q(t+1)-q(t-1))/(2*TA)*q(t)
    for ii = 2:N-1
        buf = quatmultiply(2*(x(4:7,ii+1)'-x(4:7,ii-1)')/(2*TA),quatinv(x(4:7,ii)'));
        dx(4:6,ii) = buf(2:4)';
    end
    % Boundaries: set to neighbouring value
    dx(:,1) = dx(:,2); dx(:,N) = dx(:,N-1);
end

function ddx = get_ddx(dx,TA)
% Compute the linear and angular acceleration from the velocity.
    N = size(dx,2);
    ddx = zeros(6,N);
    % Get linear velocity by computing v(t) = (x(t+1)-x(t-1))/(2*TA)
    for ii = 2:N-1
        for kk = 1:6
            ddx(kk,ii) = (dx(kk,ii+1)-dx(kk,ii-1))/(2*TA);
        end
    end
    % Boundaries: set to neighbouring value
    ddx(:,1) = ddx(:,2); ddx(:,N) = ddx(:,N-1);
end

function [dx_filt, ddx_filt] = derive_and_mean_filt(x, k, TA)
% Applies derivation and moving average filtering to calculate dx and ddx
% from x. To specify multiple filter window sizes, k can be a vector.     
    n_filt = length(k);
    n = size(x,2);
    dx_filt = zeros(6*n_filt, n);
    ddx_filt = zeros(6*n_filt, n);
    dx_filt_buf = get_dx(x, TA);
    for ii = 1:n_filt
        dx_filt((1:6)+6*(ii-1), :) = movmean(dx_filt_buf,k(ii),2);
        ddx_filt_buf = get_ddx(dx_filt((1:6)+6*(ii-1), :), TA);
        ddx_filt((1:6)+6*(ii-1), :) = movmean(ddx_filt_buf, k(ii),2);
    end
end


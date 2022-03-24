%%This function produces the known input and output terms of the
%%translational model from the decoupled estimator. !NOT ADAPTED FOR
%%INHOMOGENEOUS ROBOT TEAMS!
function [Phi,target] = model_trans(n,g,m,d,k,ddx_i,dx_i,x_i,qi_init,ddx_i_des,dx_i_des,x_i_des,ddx_j_des,dx_j_des,x_j_des)
% local states
ddp_i = ddx_i(1:3)';
dw_i = ddx_i(4:6)';
dp_i = dx_i(1:3)';
w_i = dx_i(4:6)';
p_i = x_i(1:3)';
q_i = x_i(4:7);

% all desired states
% ddp_i_des = ddx_i_des(1:3)';
% dp_i_des = dx_i_des(1:3)';
ddp_i_des = zeros(3,1);
dp_i_des = zeros(3,1);
p_i_des = x_i_des(1:3)';

% ddp_j_des = ddx_j_des(1:3)';
% dp_j_des = dx_j_des(1:3)';
ddp_j_des = zeros(3,1);
dp_j_des = zeros(3,1);
p_j_des = x_j_des(1:3)';

% skew-symmetric matrices and caligraphic S
T = findCalS(dw_i,w_i);
S_wi = findSkew(w_i);

% rotation matrix
q_delta = quatmultiply(q_i,quatinv(qi_init));
w_R_o = quat2rotm(q_delta);


target = -ddp_i*(2*m) + m*sum([ddp_i_des,ddp_j_des],2) - dp_i*2*d +...
    d*sum([dp_i_des,dp_j_des],2) - p_i*2*k + k*sum([p_i_des,p_j_des],2) ;

% Input side calculation, see eq.(15)
Phi = [T*w_R_o*makeWideMat([m,m])+S_wi*w_R_o*makeWideMat([d,d])+w_R_o*makeWideMat([k,k]),-T*w_R_o,ddp_i-g];
Phi(:,3*(n-1)+1:3*n) = [];

end
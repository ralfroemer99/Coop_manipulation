%%This function produces the known input and output terms of the rotational model from the decoupled estimator. 
% See thesis, section 4.3.2
% !NOT ADAPTED FOR INHOMOGENEOUS ROBOT TEAMS!
% function [Phi,target] = model_rot(ri_hat,rj_hat,m,d,k,j_eef,delta,kappa,ddx_i,dx_i,x_i,qi_init,qj_init,ddx_i_des,dx_i_des,x_i_des,ddx_j_des,dx_j_des,x_j_des)
function [Phi,target,target2,imp_l,imp_r] = model_rot(ri_hat,rj_hat,m,d,k,j_eef,delta,kappa,ddx_i,dx_i,x_i,qi_init,qj_init,ddx_i_des,dx_i_des,x_i_des,ddx_j_des,dx_j_des,x_j_des,ddx_j,dx_j,x_j,h_i,h_j)
% States of the local EEF.
p_i = x_i(1:3)';
q_i = x_i(4:7);
dp_i = dx_i(1:3)';
w_i = dx_i(4:6)';
ddp_i = ddx_i(1:3)';
dw_i = ddx_i(4:6)';

% Rotation matrix of object frame relative to the world frame.
q_delta = quatmultiply(q_i,quatinv(qi_init));
R_o = quat2rotm(q_delta);

% transform r from object to world frame
r_i = R_o*ri_hat;
r_j = R_o*rj_hat;

% Transformation to other Eef.
% p_j = p_i + (r_j-r_i);
% dp_j = dp_i + findSkew(w_i)*(r_j-r_i);
% ddp_j = ddp_i + (findCalS(dw_i,w_i))*(r_j-r_i);

% REMOVE: %%%%%%%%%%
p_j = x_j(1:3)';
dp_j = dx_j(1:3)';
ddp_j = ddx_j(1:3)';
%%%%%%%%%%%%%%%%%%%%

% Skew Matrices
S_ri = findSkew(r_i);
S_rj = findSkew(r_j);

% Desired states
% ddp_i_des = ddx_i_des(1:3)';
% dw_i_des = ddx_i_des(4:6)';
% ddp_j_des = ddx_j_des(1:3)';
% 
% dp_i_des = dx_i_des(1:3)';
% w_i_des = dx_i_des(4:6)';
% dp_j_des = dx_j_des(1:3)';

% REMOVE: %%%%%%%%%%%%%%%%%%%%%%%%
ddp_i_des = zeros(3,1);
dw_i_des = zeros(3,1);
ddp_j_des = zeros(3,1);

dp_i_des = zeros(3,1);
w_i_des = zeros(3,1);
dp_j_des = zeros(3,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_i_des = x_i_des(1:3)';
p_j_des = x_j_des(1:3)';

q_i_des = x_i_des(4:7);
q_j_des = x_j_des(4:7);

% Eq. (4.69) - (4.71) --> correctly implemented
R11 = R_o(1,1);
R12 = R_o(1,2);
R13 = R_o(1,3);
R21 = R_o(2,1);
R22 = R_o(2,2);
R23 = R_o(2,3);
R31 = R_o(3,1);
R32 = R_o(3,2);
R33 = R_o(3,3);
V =  [R11^2,   2*R11*R12,       2*R11*R13,       R12^2,   2*R12*R13,       R13^2;
        R11*R21, R11*R22+R21*R12, R13*R21+R11*R23, R12*R22, R12*R23+R13*R22, R13*R23;
        R11*R31, R12*R31+R11*R32, R11*R33+R13*R31, R12*R32, R12*R33+R13*R32, R13*R33;
        R21^2,   2*R21*R22,       2*R21*R23,       R22^2,   2*R22*R23,       R23^2;
        R21*R31, R21*R32+R22*R31, R21*R33+R23*R31, R22*R32, R22*R33+R32*R23, R23*R33;
        R31^2,   2*R31*R32,       2*R31*R33,       R32^2,   2*R32*R33,       R33^2];

% Eq. (4.21) in Master thesis --> correct
W = [dw_i(1)        , dw_i(2)-w_i(1)*w_i(3)  , dw_i(3)+w_i(1)*w_i(2)  , -w_i(2)*w_i(3)  , w_i(2)^2-w_i(3)^2      , w_i(2)*w_i(3);
     w_i(1)*w_i(3)  , dw_i(1)+w_i(2)*w_i(3)  , w_i(3)^2-w_i(1)^2      , dw_i(2)         , dw_i(3)-w_i(1)*w_i(2)  , -w_i(1)*w_i(3);
    -w_i(1)*w_i(2)  , w_i(1)^2-w_i(2)^2      , dw_i(1)-w_i(2)*w_i(3)  , w_i(1)*w_i(2)   , dw_i(2)+w_i(1)*w_i(3)  , dw_i(3)        ];
Phi = W*V;

% q_j = quatmultiply(q_i,quatinv(quatmultiply(qi_init,quatinv(qj_init)))); --> WRONG
% q_j = quatmultiply(quatmultiply(q_i, quatinv(qi_init)), qj_init);

% REMOVE: %%%%%
q_j = x_j(4:7);
%%%%%%%%%%%%%%%

delta_q_1 = quatmultiply(q_i,quatinv(q_i_des));
delta_q_2 = quatmultiply(q_j,quatinv(q_j_des));

delta_p_i = p_i - p_i_des;
delta_dp_i = dp_i - dp_i_des;
delta_ddp_i = ddp_i - ddp_i_des;

delta_p_j = p_j - p_j_des;
delta_dp_j = dp_j - dp_j_des;
delta_ddp_j = ddp_j - ddp_j_des;

% See eq. (42) in paper
target = - S_ri * (m*delta_ddp_i + d*delta_dp_i + k*delta_p_i) - S_rj * (m*delta_ddp_j + d*delta_dp_j + k*delta_p_j) ...
         - 2*j_eef*(dw_i - dw_i_des) - 2*delta*(w_i - w_i_des) - 2*kappa*(delta_q_1(1)*delta_q_1(2:4)' + delta_q_2(1)*delta_q_2(2:4)');

% REMOVE: %%%%%%
f_i = h_i(1:3)';
tau_i = h_i(4:6)';
f_j = h_j(1:3)';
tau_j = h_j(4:6)';

target2 = -S_ri*f_i - tau_i - S_rj*f_j - tau_j;         
%%%%%%%%%%%%%%%%

imp_l = [(m*delta_ddp_i + d*delta_dp_i + k*delta_p_i); j_eef*(dw_i - dw_i_des) + delta*(w_i - w_i_des) + 2*kappa*(delta_q_1(1)*delta_q_1(2:4)')];
imp_r = [f_i; tau_i];

end

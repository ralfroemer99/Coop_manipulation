%%Function that resolves the consensus states that are used in order to
%%calculate the consensus of the CoM and the respective covariance. 
%%See eq.(34)
function [mu, Sigma] = ResolveChi(xi)
    mu = zeros(7,1);
    Sigma = zeros(7,1);
    for ii = 1:7
        mu(ii) = xi(ii)/xi(7+ii);
        Sigma(ii) = 1/xi(7+ii);
    end
%     r_hat = zeros(3,1);
%     r_hat(1) = xi_r(1)*1/((xi_r(4)));
%     r_hat(2) = xi_r(2)*1/((xi_r(5)));
%     r_hat(3) = xi_r(3)*1/((xi_r(6)));
%     Sigma_r_hat = zeros(3,1);
%     Sigma_r_hat(1) = 1/(xi_r(4));
%     Sigma_r_hat(2) = 1/(xi_r(5));
%     Sigma_r_hat(3) = 1/(xi_r(6));
end
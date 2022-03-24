function [ri_hat, rj_hat, mo_hat, ri_tilde, rj_tilde, mo_tilde] = ratioDistribution(mu, Sigma)
if size(Sigma, 1) == size(Sigma, 2)
    Sigma = diag(Sigma);
end

    rji_hat = mu(1:3);
    mo0ri_hat = mu(4:6);
    mo_hat = mu(7);
    
    rji_tilde = Sigma(1:3);
    mo0ri_tilde = Sigma(4:6);
    mo_tilde = Sigma(7);
    
    ri_hat = mo0ri_hat./mo_hat;     % eq. (25)
    rj_hat = ri_hat + rji_hat;      % eq. (27)
    
    ri_tilde = zeros(3,1);
    rj_tilde = zeros(3,1);
    for ii = 1:3
        ri_tilde(ii) = ri_hat(ii)^2*(mo0ri_tilde(ii)/(mo0ri_hat(ii)^2) + mo_tilde/(mo_hat^2));      % eq. (26)
        rj_tilde(ii) = ri_tilde(ii) + rji_tilde(ii);                                                % eq. (28)
    end
    
    ri_tilde = diag(ri_tilde);
    rj_tilde = diag(rj_tilde);
    
%     mu_hat = [ri_hat', rj_hat', mo_hat]';
%     Sigma_hat = [ri_tilde', rj_tilde', mo_tilde]';
end



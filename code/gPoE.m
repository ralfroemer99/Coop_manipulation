%%This function calculates the generalized product of experts for a
%%variable amount of inputs (input request very basic, should be adapted).
% See eq.(23),(24)
function [mu_gpoe, Sigma_gpoe] = gPoE(mu, Sigma)
    if size(mu,1) ~= size(Sigma,1)
        error('Check dimensions of mu and Sigma!');
    end
    m = size(mu,2);
    
    switch ndims(Sigma)
        % Scalar random variable, i.e. size(mu) = size(Sigma) = (1,m)
        case 2
            Sigma_gpoe = 0;
            for ii = 1:m
                Sigma_gpoe = Sigma_gpoe + 1/Sigma(ii);
            end
            Sigma_gpoe = m*Sigma_gpoe^(-1);
            
            mu_gpoe = 0;
            for ii = 1:m
                mu_gpoe = mu_gpoe + Sigma(ii)^(-1)*mu(ii);
            end
            mu_gpoe = 1/m*Sigma_gpoe*mu_gpoe;
        % n-dimensional random variable, i.e. size(mu) = (n,m), size(Sigma) = (n,n,m)
        case 3
            Sigma_gpoe = zeros(size(Sigma,1), size(Sigma,2));
            for ii = 1:m
                Sigma_gpoe = Sigma_gpoe + Sigma(:,:,ii)^(-1);
            end
            Sigma_gpoe = m*Sigma_gpoe^(-1);
    
            mu_gpoe = zeros(size(mu,1), 1);
            for ii = 1:m
                mu_gpoe = mu_gpoe + Sigma(:,:,ii)^(-1)*mu(:,ii);
            end
    
            mu_gpoe = 1/m*Sigma_gpoe*mu_gpoe;
    end
end
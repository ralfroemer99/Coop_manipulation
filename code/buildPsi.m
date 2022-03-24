function psi = buildPsi(mu, Sigma)
    n = length(mu);
    if size(Sigma, 1) == size(Sigma, 2)
        Sigma = diag(Sigma);
    end
    psi = zeros(2*n, 1);
    for ii = 1:n
        psi(ii) = mu(ii)/Sigma(ii);
        psi(ii+n) = 1/Sigma(ii);
    end
end


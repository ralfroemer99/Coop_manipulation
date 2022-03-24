%%Function that resolves the consensus states that are used in order to
%%calculate the consensus of the inertia and the respective covariance
function [jo_hat, Sigma_jo] = ResolveChiJo(chi_jo)
    Sigma_vec = chi_jo(7:42);
    Sigma_jo = 1./reshape(Sigma_vec,[6,6]);
    jo_hat = Sigma_jo*chi_jo(1:6);
    Sigma_jo = diag(Sigma_jo);
end

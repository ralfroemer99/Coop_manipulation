function [xi_i_new] = dynAvgCons(i, xi_i, xi_j, psi_i, psi_i_old, A)
    if i == 1
        j = 2;
    else
        j = 1;
    end
    xi_i_new = xi_i + A(i,j)*(xi_j - xi_i) + psi_i - psi_i_old;
end

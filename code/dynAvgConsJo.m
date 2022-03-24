%%Function that calculates the states of the dyn. avg. cons. for the
%%inertia
%%See eq.(31)-(33)
function chi_i_new  = dynAvgConsJo(i,xi_i,xi_j,jo_hat,jo_hat_old,Sigma_jo_hat,Sigma_jo_hat_old,A)
    if i == 1
        j = 2;
    else
        j = 1;
    end

    a = (Sigma_jo_hat)^(-1)*jo_hat;
    b = (Sigma_jo_hat_old)^(-1)*jo_hat_old;
    c = 1./(reshape(Sigma_jo_hat,[36,1]));
    d = 1./(reshape(Sigma_jo_hat_old,[36,1]));
    
    chi_i_new = (xi_i - A(i,j)*(xi_i-xi_j)  + [a;c] - [b;d]);    
end



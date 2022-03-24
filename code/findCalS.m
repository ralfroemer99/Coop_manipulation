%%Calculate the specific matrix for transforming the constraints of
%%accelerations in the interaction dynamics model 
function S = findCalS(dw,w)
    S = findSkew(dw) + (findSkew(w))^2;
end


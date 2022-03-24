%%Calculates the skew symmetric matrix S(.) that calculates the
%%3-dimensional cross-product
function S = findSkew(r)
    S = [ 0,    -r(3),  r(2);
          r(3),  0,    -r(1);
         -r(2),  r(1),  0];
end
%%This function makes the wide matrices as needed for the interaction
%%dynamics model (setting 3x3 matrices side by side)
function M = makeWideMat(vec)
    n = length(vec);
    M = zeros(3,3*n);
    for i = 1:n
       M(1:3,(3*(i-1)+1):3*i) = vec(i)*eye(3);
    end
end


function j_guess = guessInertia(m_o, n)
% Estimate the object inertia tensor expressed in the object frame by approximating it with a cuboid of homogeneous density.
% m_o:  Object mass without additional weights
% n:    Number of additional weights added to the object. Each has a mass of 1.043 kg.

% Dimensions
x = 0.3;
y = 0.6;
z = 0.3;

% Additional weight
m_add = 1.043;

% Calculate inertia tensor of the object without additional weights
jxx = m_o/12*(y^2 + z^2);
jxy = 0;
jxz = 0;
jyy = m_o/12*(x^2 + z^2);
jyz = 0;
jzz = m_o/12*(x^2 + y^2);

% Add weights
% Distance to x-axis
dx1 = sqrt((0.3-0.083)^2 + 0.045^2);
dx2 = sqrt((0.3-0.083)^2 + 0.09^2);

% Distance to y-axis
dy = 0.09;

% Distance to z-axis
dz1 = sqrt((0.3-0.083)^2 + 0.07942^2);
dz2 = 0.3-0.083;


jxx = jxx + 4*m_add*dx1^2 + 2*m_add*dx2^2;
jyy = jyy + 6*m_add*dy^2;
jzz = jzz + 4*m_add*dz1^2 + 2*m_add*dz2^2;

j_guess = [jxx; jxy; jxz; jyy; jyz; jzz];

end
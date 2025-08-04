function rtn = quat2D(quat)
% Initiate
a = quat(1);
b = quat(2);
c = quat(3);
d = quat(4);

% Define rotation matrix
rtn(1,1) = a*a + b*b - c*c - d*d;
rtn(1,2) = 2*b*c - 2*a*d;
rtn(2,1) = 2*b*c + 2*a*d;
rtn(2,2) = a*a - b*b + c*c - d*d;
end
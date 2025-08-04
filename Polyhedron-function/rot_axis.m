% Thi VO
% Columbia University
% Point rotator about arbitrary axis

function rot_mat = rot_axis(n_vect,theta)
% Input
% n_vect: vector of axis
% theta: angle to rotate
% point: point to be rotated

% Output:
% point_rot: rotated point

% Parameters defining axis
l = n_vect(1);
m = n_vect(2);
n = n_vect(3);

% Creating transformation matrix
row_1 = [(l*l*(1-cos(theta))+cos(theta)) (m*l*(1-cos(theta))-n*sin(theta)) (n*l*(1-cos(theta))+m*sin(theta))];
row_2 = [(l*m*(1-cos(theta))+n*sin(theta)) (m*m*(1-cos(theta))+cos(theta)) (n*m*(1-cos(theta))-l*sin(theta))];
row_3 = [(l*n*(1-cos(theta))-m*sin(theta)) (m*n*(1-cos(theta))+l*sin(theta)) (n*n*(1-cos(theta))+cos(theta))];

% Combining transformation matrix
rot_mat = [row_1; row_2; row_3];

end
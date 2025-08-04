% Thi Vo
% Columbia University
% Rotation Matrix - Euler's Angles

function [point_rot] = rot_euler(psi,theta,phi)
% Inputs
% point: point to be rotated
% psi: angle between x and N axis
% theta: angle between z and Z axis
% phi: angle between N and X axis

% Output:
% point_rot: rotated_point

% Defining rotation matrix

rot_matx = [cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi) sin(psi)*cos(theta)*cos(phi)+cos(psi)*sin(phi) -sin(theta)*cos(phi);
    -cos(psi)*cos(theta)*sin(phi)-sin(psi)*cos(phi) -sin(psi)*cos(theta)*sin(phi)+cos(psi)*cos(phi) sin(theta)*sin(phi);
    cos(psi)*sin(theta) sin(psi)*sin(theta) cos(theta)];

% Performing roation
point_rot = rot_matx;

function eul_angle = rtx_to_euler()

sy = sqrt(R(1,1)^2 + R(2,1)^2);

% Two cases
if sy > 1E-6
    z = atan2(R(3,2) , R(3,3));
    y = atan2(-R(3,1), sy);
    x = atan2(R(2,1), R(1,1));
else
    z = atan2(-R(2,3), R(2,2));
    y = atan2(-R(3,1), sy);
    x = 0;
end

eul_angle = [-x, -y, -z];
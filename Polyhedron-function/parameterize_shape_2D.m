function [pos_theo,theta_space_poly,Fx,Fy] = parameterize_shape_2D(pts_poly,grid_step)
% Convert to spherical
theta_poly = atan2(pts_poly(:,2),pts_poly(:,1));
% Make space
d_grid = grid_step;
shift_value = 0;
theta_space_poly = -pi+shift_value:d_grid:pi-shift_value;
% Map surface
method = 'linear';
extend = 'extrap';
Fx = interp1(theta_poly,pts_poly(:,1),theta_space_poly',method,extend);
Fy = interp1(theta_poly,pts_poly(:,2),theta_space_poly',method,extend);
% Store
pos_theo = [Fx Fy];
end

clear; close all; clc;

warning off;
addpath('C:\Users\lanth\Documents\codes\polyhedron_function')

%%%%%%%%%%%%%%%%%
%%% Load data %%%
%%%%%%%%%%%%%%%%%

% Scale factor
corona_factor = 0.45;

% Save name
file_name_data = 'rhombic_sample_corona.mat'

% Load data
load(file_name_data)

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%

%%% Define colors %%%
purple = [107.0/255, 62.0/255, 154.0/255];
blue = [0, 64.0/255, 128.0/255];
red = [180.0/255, 7.0/255, 5.0/255];
green = [31 164 10]/255;
cyan = [124 207 247]/255;
colors = {red,blue,purple,green};

% Paper colors
core_color = [206 186 37]/255;
core_color = [211 207 209]/255;
corona_color = [171 202 247]/255;
iodine_color = [108 0 214]/255;
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Interpolate and make mesh %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define colormap
% Blue
color0 = [0 25 51];
color1 = 1.15*[0 76 153];
color2 = 1.25*[0 102 204];
color3 = [0 128 255];
color4 = [51 153 255];
color5 = [102 179 255];
color6 = [153 204 255];
color7 = [192 192 192];
color8 = [224 224 224];
color_mtx = [color8; color7; color6; color5; color4; color3; color2; color1; color0]/255;
% Generate grid
indx_plot = 1:1:n_graft;
n_color = 1.1*length(indx_plot);
n_color_type = length(color_mtx(:,1));
n_color_grid = ceil(n_color/n_color_type);
mymap_corona = [];
for i = 1:n_color_type-1
    tmp_ith = color_mtx(i,:);
    tmp_jth = color_mtx(i+1,:);
    tmp_mtx = [linspace(tmp_ith(1),tmp_jth(1),n_color_grid)', linspace(tmp_ith(2),tmp_jth(2),n_color_grid)', linspace(tmp_ith(3),tmp_jth(3),n_color_grid)'];
    mymap_corona = [mymap_corona; tmp_mtx];
end

% Purple
color0 = [25 0 51];
color1 = [51 0 102];
color2 = [76 0 153];
color3 = [102 0 204];
color4 = [127 0 255];
color5 = [153 51 255];
color6 = [178 105 255];
color7 = [192 192 192];
color8 = [224 224 224];
color_mtx = [color7; color6; color5; color4; color3; color2; color1]/255;
% Generate grid
indx_plot = 1:1:n_graft;
n_color = 1.1*length(indx_plot);
n_color_type = length(color_mtx(:,1));
n_color_grid = ceil(n_color/n_color_type);
mymap_iodine = [];
for i = 1:n_color_type-1
    tmp_ith = color_mtx(i,:);
    tmp_jth = color_mtx(i+1,:);
    tmp_mtx = [linspace(tmp_ith(1),tmp_jth(1),n_color_grid)', linspace(tmp_ith(2),tmp_jth(2),n_color_grid)', linspace(tmp_ith(3),tmp_jth(3),n_color_grid)'];
    mymap_iodine = [mymap_iodine; tmp_mtx];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generate full patch shape %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute norm
dr_core = sqrt(sum(pts_core.^2,2));
dr_corona = sqrt(sum(pts_corona.^2,2));

% Compute displacement
dr_diff = dr_corona - dr_core;
indx_diff = find(dr_diff > -1E-2);

% Omega
omega_core = dr_core;

% Generate sphere
pts_sph = points_on_sphere(50,1);

% Loop and compute net growth
pts_corona_full = [];
for i = 1:length(indx_diff)
    Rgraft_tmp = (R_graft(indx_diff(i)) - 0*ro)/ro;
    Rgraft_growth = 1;
    xi_r = 0;
    Rgraft_pos = 0;
    while Rgraft_growth < Rgraft_tmp
        % Compute blob size
        xi_tmp = Rgraft_growth/sqrt(n_graft)*omega_core(indx_diff(i)).^(3/2);
        % Update
        Rgraft_growth = Rgraft_growth + xi_tmp;
        % Store
        Rgraft_pos(end+1) = Rgraft_growth;
        xi_r(end+1) = 3*xi_tmp;
    end
    % Update
    Rgraft_pos = Rgraft_pos(2:end);
    xi_r = xi_r(2:end);
    % Add to mesh
    for ii = 1:length(xi_r)
        % Draw blob
        pts_xi_tmp = xi_r(ii)*pts_sph;
        % Direction
        nv_growth = pts_core(indx_diff(i),:);
        nv_growth = nv_growth/norm(nv_growth);
        % Place
        pts_xi_tmp = bsxfun(@plus,pts_xi_tmp,pts_core(indx_diff(i),:)+nv_growth*Rgraft_pos(ii));        
        % Update
        pts_corona_full = [pts_corona_full; pts_xi_tmp];
    end
end
shrink_factor = 0.5;
pts_corona_full = shrink_factor*pts_corona_full;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
%%% Plotting %%%
%%%%%%%%%%%%%%%%

figure(1)
hold on;
core_plot = plot(alphaShape(pts_verts,1E10),'EdgeColor','none','FaceColor',core_color);
pts_corona_full = corona_factor*pts_corona_full;
shp_corona = alphaShape(1.05*pts_corona_full,0.3);
corona_plot = plot(shp_corona,'EdgeColor','none','FaceAlpha',0.675','FaceColor',cyan);
for j = 1:length(faces)
    tmp = pts_verts;
    tmp_faces = faces{j};
    tmp_faces = [tmp_faces tmp_faces(1)];
    plot3(tmp(tmp_faces,1),tmp(tmp_faces,2),tmp(tmp_faces,3),'k','Linewidth',3.5)
end
view_angle = [-90 22];
camlight
lighting gouraud
set(core_plot,'facelighting','flat')
set(corona_plot,'facelighting','gouraud')
axis equal
axis off
xlabel('x')
ylabel('y')
zlabel('z')
camlight(view_angle(1)+45,view_angle(2)-30)
view(view_angle)
camproj('perspective')

%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
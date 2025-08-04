clear; close all; clc;

warning off;
addpath('C:\Users\lanth\Documents\codes\polyhedron_function')

%%% Colors %%%
purple = [107.0/255, 62.0/255, 154.0/255];
blue = [0, 64.0/255, 128.0/255];
red = [180.0/255, 7.0/255, 5.0/255];
green = [0.0/255, 112/255, 0.0/255];
yellow = [241, 194, 50]/255;
colors = {red, blue, purple, green, yellow};
%%%%%%%%%%%%%%

%%% Load mesh %%%
% Define shape
verts_name = 'RhombicDodecahedron';
mesh_name = strcat(lower(verts_name),'_mesh.mat');
% Load
mesh_data = load(mesh_name);
omega_original = mesh_data.omega;
% Current directory
current_path = pwd;
%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PARAMETERS TO CHANGE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Chemical potential (controls concentration [I])
mu = -0.175;

% Grafting density
fo = 0.1;

% Chain-chain attraction
chi = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Determine grafting ratio %%%
% Energy
eps_A = 1;
eps_B = 1;
eps_C = -1;
% Surface packing
rho_ref = 0.9;
rho_A = 0.79/rho_ref;
rho_B = 0.56/rho_ref;
rho_C = 0.9/rho_ref;
% Power Scaling
ppower = 3/3;
rho_A = rho_A.^ppower;
rho_B = rho_B.^ppower;
rho_C = rho_C.^ppower;
% Temperature
kT = 0.5;
% Flag for removing certain facets
if strcmp(verts_name,'Octahedron') == 1
    flag_A = 1;
    flag_B = 0;
    flag_C = 1;
elseif strcmp(verts_name,'RhombicDodecahedron') == 1
    flag_A = 0;
    flag_B = 1;
    flag_C = 1;
elseif strcmp(verts_name,'Cuboctahedron') == 1
    flag_A = 1;
    flag_B = 0;
    flag_C = 1;
elseif strcmp(verts_name,'trunc_tetra') == 1
    flag_A = 1;
    flag_B = 0;
    flag_C = 1;
elseif strcmp(verts_name,'Dipyramid5') == 1
    flag_A = 1;
    flag_B = 0;
    flag_C = 1;
end
% Working terms
D_A = exp( flag_A*(1/kT).*(eps_C.*rho_A - mu) );
D_B = exp( flag_B*(1/kT).*(eps_C.*rho_B - mu) );
D_C = exp( flag_C*(1/kT).*(eps_C.*rho_C - mu) );
% A Type Fraction
phi_A1 = 1 - D_A./(1+D_B+D_B.*D_A);
phi_A2 = D_C + D_A./(1+D_A).*(1 - D_A./(1+D_B+D_B.*D_A));
phi_A = D_C./(1+D_A).*phi_A1./phi_A2;
% B Type Fraction
phi_B1 = D_C.*D_A./(1+D_B+D_B.*D_A);
phi_B2 = D_C + D_A./(1+D_A).*(1 - D_A./(1+D_B+D_B.*D_A));
phi_B = phi_B1./phi_B2;
% C TYpe Fraction
phi_C1 = D_A./(1+D_A).*(1-D_A./(1+D_B+D_B.*D_A));
phi_C = phi_C1./(D_C+phi_C1);
% Initial
phi_Tinit = phi_A + phi_B + phi_C;
% Check flags
if flag_A == 0
    phi_A = 0*phi_A;
elseif flag_B == 0
    phi_B = 0*phi_B;
elseif flag_C == 0
    phi_C = 0*phi_C;
end  
% Combine
phi_T = phi_A + phi_B + phi_C;
% Define relative probability
if flag_A == 0
    phi_A = 0*phi_A;
    % Scale
    phi_Brel = phi_B/(phi_B+phi_C);
    phi_Crel = phi_C/(phi_B+phi_C);
    % Define working probability
    phi_0 = phi_Crel;
    phi_1 = phi_Brel;
elseif flag_B == 0
    phi_B = 0*phi_B;
    % Scale
    phi_Arel = phi_A/(phi_A+phi_C);
    phi_Crel = phi_C/(phi_A+phi_C);
    % Define working probability
    phi_0 = phi_Crel;
    phi_1 = phi_Arel;
elseif flag_C == 0
    phi_C = 0*phi_C;
    % Scale
    phi_Arel = phi_A/(phi_A+phi_B);
    phi_Brel = phi_B/(phi_A+phi_B);
    % Define working probability
    phi_0 = phi_Arel;
    phi_1 = phi_Brel;
end  
[phi_A phi_B phi_C phi_T phi_0 phi_1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Flag for removing certain facets %%%
if strcmp(verts_name,'Octahedron') == 1
    faces = {[1 2 4],[1 5 4],[5 6 4],[2 6 4],...
        [1 2 3],[1 5 3],[5 6 3],[2 6 3]};
elseif strcmp(verts_name,'RhombicDodecahedron') == 1
    faces = {[6 8 13 11],[6 8 3 1],[3 8 10 4],[10 8 13 14],...
        [11 13 14 12],[10 4 9 14],[4 3 1 2],[1 6 11 5],...
        [5 11 12 7],[12 14 9 7],[7 9 4 2],[7 2 1 5]};
elseif strcmp(verts_name,'Cuboctahedron') == 1
    faces_alt = {[3 6 9],[2 6 8],[9 11 12],[8 12 10],[11 5 7],[10 7 4],[5 1 3],[4 1 2]};
    faces = {[3 5 11 9],[11 12 10 7],[1 5 7 4],[1 3 6 2],[6 9 12 8],[2 8 10 4]};       
elseif strcmp(verts_name,'Dipyramid5') == 1
    faces = {[2 6 4],[2 4 3],[2 3 5],[2 5 7],[2 7 6],...
        [1 6 4],[1 4 3],[1 3 5],[1 5 7],[1 7 6]};    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Assign facet to shape %%%
% Cutoff
phi_T_cutoff = 0.95;
% Define terms
pts_verts = mesh_data.pts_verts;
pts_mesh = mesh_data.pts_mesh;
if roundn(phi_T,-2) < phi_T_cutoff
    if (strcmp(verts_name,'Cuboctahedron') == 1) || (strcmp(verts_name,'trunc_tetra') == 1)
        
        %%% Define reference %%%
        face_indx_ref = zeros(length(pts_mesh(:,1)),1);
        indx_full = [];
        for i = 1:length(faces)
            % Define direction
            pts_tmp = pts_verts(faces{i},:);
            nv_tmp = mean(pts_tmp);
            nv_tmp = nv_tmp/norm(nv_tmp);
            % Define rotation matrix
            rtx_tmp = rot_az(nv_tmp,1);
            % Rotate mesh
            pts_mesh_rot = pts_mesh*rtx_tmp;
            % Rotate face
            pts_face_rot = pts_tmp*rtx_tmp;
            % Generate polyhedron
            shift_value = 0.05;
            pts_region_tmp1 = pts_face_rot;
            pts_region_tmp1(:,3) = pts_region_tmp1(:,3) - shift_value;
            pts_region_tmp2 = pts_face_rot;
            pts_region_tmp2(:,3) = pts_region_tmp2(:,3) + shift_value;
            pts_region_tmp = [pts_region_tmp1; pts_region_tmp2];
            pts_region_shape = alphaShape(pts_region_tmp,1E10);
            % Check in shape
            indx_in = inShape(pts_region_shape,pts_mesh_rot);
            indx_in = find(indx_in == 1);
            indx_full = [indx_full; indx_in];
        end
        % Set face
        indx_full = unique(indx_full);
        face_indx_ref(indx_full) = 1;
        % Compute max fraction
        indx_tmp = find(face_indx_ref == 1);
        max_frac_coverage = length(indx_tmp)/length(face_indx_ref);
        %%%%%%%%%%%%%%%%%%%%%%%% 
        
        % Initialize
        face_indx = zeros(length(pts_mesh(:,1)),1);
        
        % Compute kernel
        kernel_tmp = sqrt(sum(pts_mesh.^2,2));
        kernel_tmp = kernel_tmp./min(kernel_tmp(:));
        
        % Check if grow or shrink
        if max_frac_coverage > phi_T
            %%%%%%%%%%%%%%%%%%%%%%%
            %%% Shrink coverage %%%
            %%%%%%%%%%%%%%%%%%%%%%%
            % Scale verts
            verts_scale_factor = 0.99;
            verts_scale = 1/verts_scale_factor;
            % Set flag
            flag_scale = 1;
            while flag_scale == 1
                % Reset
                indx_full = [];
                % Loop
                for i = 1:length(faces)
                    % Define direction
                    pts_tmp = pts_verts(faces{i},:);
                    nv_tmp = mean(pts_tmp);
                    nv_tmp = nv_tmp/norm(nv_tmp);
                    % Define rotation matrix
                    rtx_tmp = rot_az(nv_tmp,1);
                    % Rotate mesh
                    pts_mesh_rot = pts_mesh*rtx_tmp;
                    % Rotate face
                    pts_face_rot = pts_tmp*rtx_tmp;
                    % Scale
                    verts_scale = verts_scale*verts_scale_factor;
                    pts_face_rot(:,1:2) = verts_scale*pts_face_rot(:,1:2);
                    % Generate polyhedron
                    shift_value = 0.05;
                    pts_region_tmp1 = pts_face_rot;
                    pts_region_tmp1(:,3) = pts_region_tmp1(:,3) - shift_value;
                    pts_region_tmp2 = pts_face_rot;
                    pts_region_tmp2(:,3) = pts_region_tmp2(:,3) + shift_value;
                    pts_region_tmp = [pts_region_tmp1; pts_region_tmp2];
                    % Create alphaShape
                    pts_region_shape = alphaShape(1.01*pts_region_tmp,1E10);
                    % Check in shape
                    indx_in = inShape(pts_region_shape,pts_mesh_rot);
                    indx_in = find(indx_in == 1);
                    % Store
                    indx_full = [indx_full; indx_in];
                    
                end
                % Reduce
                indx_full = unique(indx_full);
                % Compare
                frac_compare = length(indx_full)/length(kernel_tmp);
                if frac_compare < phi_T
                    flag_scale = 0;
                end
            end
            % Assign
            face_indx(indx_full) = 1;
            %%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%
        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Grow into other face region %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Scale verts
            verts_scale_factor = 0.99;
            verts_scale = 1/verts_scale_factor;
            % Set flag
            flag_scale = 1;
            while flag_scale == 1
                % Reset
                face_indx = face_indx_ref;
                indx_full = [];
                % Loop
                for i = 1:length(faces_alt)
                    % Define direction
                    pts_tmp = pts_verts(faces_alt{i},:);
                    nv_tmp = mean(pts_tmp);
                    nv_tmp = nv_tmp/norm(nv_tmp);
                    % Define rotation matrix
                    rtx_tmp = rot_az(nv_tmp,1);
                    % Rotate mesh
                    pts_mesh_rot = pts_mesh*rtx_tmp;
                    % Rotate face
                    pts_face_rot = pts_tmp*rtx_tmp;
                    % Set base
                    shift_value = 0.05;
                    pts_region_tmp1 = pts_face_rot;
                    pts_region_tmp1(:,3) = pts_region_tmp1(:,3) - shift_value;
                    pts_region_tmp2 = pts_face_rot;
                    pts_region_tmp2(:,3) = pts_region_tmp2(:,3) + shift_value;
                    pts_region_tmp_ref = [pts_region_tmp1; pts_region_tmp2];
                    % Scale
                    verts_scale = verts_scale*verts_scale_factor;
                    pts_face_rot(:,1:2) = verts_scale*pts_face_rot(:,1:2);
                    % Generate polyhedron
                    shift_value = 0.05;
                    pts_region_tmp1 = pts_face_rot;
                    pts_region_tmp1(:,3) = pts_region_tmp1(:,3) - shift_value;
                    pts_region_tmp2 = pts_face_rot;
                    pts_region_tmp2(:,3) = pts_region_tmp2(:,3) + shift_value;
                    pts_region_tmp = [pts_region_tmp1; pts_region_tmp2];
                    % Create alphaShape
                    pts_region_shape_ref = alphaShape(1.01*pts_region_tmp_ref,1E10);
                    pts_region_shape = alphaShape(1.01*pts_region_tmp,1E10);
                    % Check in shape
                    indx_in_ref = inShape(pts_region_shape_ref,pts_mesh_rot);
                    indx_in_ref = find(indx_in_ref == 1);
                    indx_in = inShape(pts_region_shape,pts_mesh_rot);
                    indx_in = find(indx_in == 1);
                    indx_in_diff = setdiff(indx_in_ref,indx_in);
                    % Store
                    indx_full = [indx_full; indx_in_diff];
                    
                end
                % Reduce
                indx_full = unique(indx_full);
                % Set extra points
                face_indx(indx_full) = 1;
                indx_full = find(face_indx == 1);
                % Compare
                frac_compare = length(indx_full)/length(kernel_tmp);
                if frac_compare > phi_T
                    flag_scale = 0;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
    elseif strcmp(verts_name,'RhombicDodecahedron') == 1
        
        % Initialize
        face_indx = ones(length(pts_mesh(:,1)),1);
        
        % Compute kernel
        kernel_tmp = sqrt(sum(pts_mesh.^2,2));
        kernel_tmp = kernel_tmp./min(kernel_tmp(:));
        
        % Find 'edges'
        kernel_cutoff = 0.9;
        indx_tmp = find(kernel_tmp > kernel_cutoff);
        
        % Tune edges
        flag_scale = 1;
        kernel_cutoff_factor = 1.001;
        while flag_scale == 1
            kernel_cutoff = kernel_cutoff*kernel_cutoff_factor;
            indx_tmp = find(kernel_tmp > kernel_cutoff);
            frac_compare = length(indx_tmp)/length(kernel_tmp);
            if frac_compare < phi_T
                flag_scale = 0;
            end
        end
        
        % Assign edges
        face_indx(indx_tmp) = 0;
        
    elseif (strcmp(verts_name,'Octahedron') == 1) || (strcmp(verts_name,'Dipyramid5') == 1)
        
        % Initialize
        face_indx = zeros(length(pts_mesh(:,1)),1);
        
        % Compute kernel
        kernel_tmp = mesh_data.omega;
        
        % Find 'edges'
        kernel_cutoff = 2;
        indx_tmp = find(kernel_tmp < kernel_cutoff);
        
        % Tune edges
        flag_scale = 1;
        kernel_cutoff_factor = 0.99;
        while flag_scale == 1
            kernel_cutoff = kernel_cutoff*kernel_cutoff_factor;
            indx_tmp = find(kernel_tmp < kernel_cutoff);
            frac_compare = length(indx_tmp)/length(kernel_tmp);
            if frac_compare < phi_T
                flag_scale = 0;
            end
        end
        
        % Assign edges
        face_indx(indx_tmp) = 1;
        
    end
    
else
    face_indx = ones(length(pts_mesh(:,1)),1);
end
if (strcmp(verts_name,'Octahedron') == 1) || (strcmp(verts_name,'Dipyramid5') == 1)
    face_indx = mesh_data.face_indx;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Assign relative distance to 111 (reference plane) %%%
% Get proper faces
% 111 plane
indx_0 = find(face_indx == 0);
pts_0 = pts_mesh(indx_0,:);
% Other plane
indx_1 = find(face_indx == 1);
pts_1 = pts_mesh(indx_1,:);
% Compute distance
if (~isempty(indx_0)) && (~isempty(indx_1))
    dr_01 = pdist2(pts_0,pts_1);
    % Get minimum
    dr_01_min = zeros(length(pts_1(:,1)),1);
    for i = 1:length(pts_1(:,1))
        dr_tmp = bsxfun(@minus,pts_0,pts_1(i,:));
        dr_tmp = sqrt(sum(dr_tmp.^2,2));
        dr_01_min(i) = min(dr_tmp);
    end
    % Define distance
    omega_I = ones(size(face_indx));
    omega_I(indx_1) = 1 + dr_01_min;
    omega_weight = omega_I.^2;
else
    omega_I = ones(size(face_indx));
    omega_weight = omega_I.^2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Attach I to surface %%%
% Define number of MC runs to average
n_mc = 100;
graft_mc = zeros(size(face_indx));
for i = 1:n_mc
% Logger  
disp([i n_mc])
% Determine number of grafts
ngraft = ceil(phi_T*length(pts_mesh(:,1)));
if ngraft > length(pts_mesh(:,1))
    ngraft = length(pts_mesh(:,1));
end
% Define logging matrix
graft_logger = zeros(size(face_indx));
% Number of types
n_type_0 = length(find(face_indx == 0));
n_type_1 = length(find(face_indx == 1));
% Note: 0 is the 111 (best plane), 1 is the other plane
ncount = 0;
while ncount < ngraft  
    % Reset update
    flag_update = 0;   
    % Generate random number
    rand_val = rand;
    % Check which facet to fill
    if rand_val < phi_0
        % 111 plane
        flag_add = 0;
    else
        % Other plane (110 or 100)
        flag_add = 1;
    end
    % Pick particle to add
    indx_sel = find( (face_indx == flag_add) & (graft_logger == 0) );
    if isempty(indx_sel) 
        % When facet is filled (flip face to add)
        if flag_add == 0
            flag_add = 1;
        elseif flag_add == 1
            flag_add = 0;
        end
        % Get new index
        indx_sel = find( (face_indx == flag_add) & (graft_logger == 0) );
        % Define weights
        weight_sel = omega_weight(indx_sel);
        % Select working index
        range_val = range(weight_sel);
        if range_val > 0
            rand_range = rand;
            cutoff_tmp = min(weight_sel) + rand_range*range_val;
            indx_tmp = find(weight_sel < cutoff_tmp);
            graft_indx1 = randi([1,length(indx_tmp)]);
            graft_indx = indx_tmp(graft_indx1);
        else
            graft_indx = randi([1,length(indx_sel)]);
        end
        % Add
        graft_logger(indx_sel(graft_indx)) = 1;
        % Update
        flag_update = 1;
    else
        % Define weights
        weight_sel = omega_weight(indx_sel);
        % Select working index
        range_val = range(weight_sel);
        if range_val > 0
            rand_range = rand;
            cutoff_tmp = min(weight_sel) + rand_range*range_val;
            indx_tmp = find(weight_sel < cutoff_tmp);
            graft_indx1 = randi([1,length(indx_tmp)]);
            graft_indx = indx_tmp(graft_indx1);
        else
            graft_indx = randi([1,length(indx_sel)]);
        end
        % Add
        graft_logger(indx_sel(graft_indx)) = 1;
        % Update
        flag_update = 1;
    end       
    % Update counter
    if flag_update == 1
        ncount = ncount + 1;
    end
end
% Store
graft_mc = graft_mc + graft_logger;
end
% Average
graft_mc = graft_mc./n_mc;
% Define grafting surface probability with [I] coverage
P_I = 1 - graft_mc;
% Scale
P_I = (P_I - min(P_I))/(max(P_I)-min(P_I));
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Chain parameters %%%
% Define chain length
N = 50;
% Particle radius
ro = 35;
% Degree of truncation
trunc_value = 0;
% Kuhn length
b = 0.97;
% Excluded volume
v = (b^3);
% Grid mesh
grid_size = 30;
% Number of monte carlo runs
n_mc = 10;
pts_cutoff = 1E6;
%%%%%%%%%%%%%%%%%%%%%%%%

%%% Graft polymer %%%
% Mean number of distance
rdist_scale = mean(pdist(pts_mesh));
% Compute area
pts_verts = mesh_data.pts_verts;
area_scale = (1E-3/mesh_data.sigma_particle)^2;
ss = alphaShape(ro*pts_verts(:,1),ro*pts_verts(:,2),ro*pts_verts(:,3),1E10);
poly_area = surfaceArea(ss);
% Determing number of grafts
fo_array = (1:numel(pts_mesh(:,1)))*(4*pi*b^2)/poly_area;
[~,indx_fo] = min(abs(fo_array-fo));
n_graft = indx_fo;
% Check number of grafting points
n_free = find(P_I > 1E-2);
n_free = length(n_free);
if n_graft > n_free
    n_graft = n_free;
end
% Initialize
graft_count = zeros(numel(pts_mesh(:,1)),1);
R_graft = zeros(length(pts_mesh(:,1)),1);
% Define kernel
kernel = (omega_original).^(1);
% kernel = ones(size(omega_original));
% Identify grafting point
[~,indx_0] = max(abs(kernel.*P_I));
% Loop
for ii = 1:n_mc
    % Grafting points matrix
    graft_indx = zeros(n_graft,1);    
    for j = 1:n_graft
        fprintf('N: %i, MC Run: %i, Chain: %i out of %i\n',N,ii,j,n_graft)
        % Flag
        flag_add = 0;
        while flag_add ~= 1  
            %%% Generate grafting probability %%%
            % Initialize
            if j == 1
                omega_graft = zeros(size(kernel));
                omega_graft(indx_0) = 1
            else
                % Select for grafts
                indx_tmp = find(graft_indx ~= 0);
                pts_grafts = pts_mesh(graft_indx(indx_tmp),:);
                Rg_grafts = Rg_chain(graft_indx(indx_tmp));
                indx_full = [];
                % Find points within grafts
                for gg = 1:length(indx_tmp)
                    dr_tmp = bsxfun(@minus,pts_mesh,pts_grafts(gg,:));
                    dr_tmp = sqrt(sum(dr_tmp.^2,2))-0*rdist_scale;
                    indx_tmp = find(dr_tmp < 2.0*(Rg_grafts(gg)-ro));
                    indx_full = [indx_full; indx_tmp];
                end
                indx_full = unique(indx_full);
                % Popular omega_graft matrix
                omega_graft = zeros(size(kernel));
                omega_graft(indx_full) = 1;
            end
            % pdf handle functions
            v_chain = @(omega,omega_graft,chi) (v*sqrt((ones(size(omega))-2*chi*omega_graft).^2)).^(5/5);
            v_graft = v_chain(kernel,omega_graft,chi);
            R_chain = @(omega,N,omega_graft,chi) ro*fo^(1/5)*(v_graft).^(1/5).*b^(2/5)*N^(3/5).*omega.^(-3/5)*(b/ro)^(3/5) - 0*ro;
            R_chain_final = @(omega,N) ro*fo^(1/5)*v^(1/5)*b^(2/5)*N^(3/5)*omega.^(-3/5)*(b/ro)^(3/5) - 0.0*ro;
            R_chain_nocore = @(omega,N,omega_graft,chi) ro*fo^(1/5)*(v_graft).^(1/5).*b^(2/5)*N^(3/5).*omega.^(-3/5)*(b/ro)^(3/5) - 0.0*ro;
            P_chain = @(omega,N,omega_graft,chi) exp( -R_chain(omega,N,omega_graft,chi).^2/(N*b^2) - v_graft.*fo*ro^2*N^2./((R_chain(omega,N,omega_graft,chi)+ro).*omega).^3);
            E_chain = @(omega,N,omega_graft,chi) (-R_chain(omega,N,omega_graft,chi).^2/(N*b^2) - v*fo*ro^2*N^2./((R_chain(omega,N,omega_graft,chi)+ro).*omega).^3);            
            % Compute chain Rg
            Rg_chain = R_chain_nocore(kernel,N,omega_graft,chi)/6;           
            % Compute grafting density
            P_graft = P_chain(kernel,N,omega_graft,chi);   
            % Rescale
            graft_filled = graft_indx(graft_indx > 0);
            indx_unique = unique(graft_filled);          
            P_graft = (P_graft-min(P_graft))/(max(P_graft)-min(P_graft));
            % Combine with [I] and rescale
            P_graft = P_graft.*P_I;
            P_graft = (P_graft-min(P_graft))/(max(P_graft)-min(P_graft));
            % Subselect
            indx_free = setdiff(1:length(P_graft),graft_indx);
            P_graft_free = P_graft(indx_free);                       
            % Grafting probability
            p_min = min(P_graft_free);
            p_max = max(P_graft_free);
            p_jth = (p_max*p_min)*rand + p_min;              
            % Grabbing index
            p_indx = find(P_graft >= p_jth);
            % Select max
            if (j ~= 1) || (length(p_indx) > 1)
                P_indx_graft = P_graft(p_indx);
                [C_indx,C_ia] = setdiff(p_indx,graft_indx);
                [~,indx_graft] = max(P_graft(C_indx));
                graft_jth = C_ia(indx_graft);
                if isempty(graft_jth)
                    graft_jth = 1;
                end
            else
                % Pick random index from selected values
                graft_jth = randi([1,length(p_indx)]);
            end
            graft_jth = randi([1,length(p_indx)]);
            % Select open point
            if isempty(find(graft_indx == p_indx(graft_jth), 1))
                graft_indx(j) = p_indx(graft_jth);
                R_graft_tmp = R_chain_final(kernel,N);
                flag_add = 1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
    % Update R
    R_sum = zeros(size(R_graft_tmp));
    R_sum(graft_indx) = R_graft_tmp(graft_indx);
    R_graft = R_graft + R_sum;   
    % Update graft counters
    for j = 1:n_graft
        graft_count(graft_indx(j)) = graft_count(graft_indx(j)) + 1;
    end
end
% Get distances
graft_count = graft_count/n_mc;
R_graft = R_chain_final(omega_original,N);
R_graft = R_graft.*graft_count;
% Temp arrays
R_tmp = zeros(length(pts_mesh(:,1)),1);
for j = 1:numel(R_tmp(:,1))
    R_tmp(j) = (R_graft(j)+ro - min(R_graft+ro));
end
% Scale distance
x_total = zeros(length(pts_mesh(:,1)),1);
y_total = zeros(length(pts_mesh(:,1)),1);
z_total = zeros(length(pts_mesh(:,1)),1);
% Generate corona
for j = 1:numel(x_total)
    n_vect = [pts_mesh(j,1) pts_mesh(j,2) pts_mesh(j,3)]/norm([pts_mesh(j,1) pts_mesh(j,2) pts_mesh(j,3)]);
    x_total(j) = n_vect(1)*R_tmp(j) + pts_mesh(j,1);
    y_total(j) = n_vect(2)*R_tmp(j) + pts_mesh(j,2);
    z_total(j) = n_vect(3)*R_tmp(j) + pts_mesh(j,3);
end
% Outputs
pts_corona = [x_total y_total z_total];
%%%%%%%%%%%%%%%%%%%%%

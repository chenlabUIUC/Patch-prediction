function [pts_core,pts_corona,pdf_core,pdf_corona,sigma_particle] = gen_corona(pts_coord,vargin)

%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%
% Particle radius
ro = vargin(1);
% Degree of truncation
trunc_value = vargin(2);
% Chain length
N = vargin(3);
% Grafting density
fo = vargin(4);
% Kuhn length
beta_kuhn = vargin(5);
% Excluded volume
v = (beta_kuhn^3);
% Grid mesh
grid_size = vargin(6)+2;
% Number of monte carlo runs
n_mc = vargin(7);
% Sigma step size
dsigma = 0.05*ro;
% Cutoff points
pts_cutoff = vargin(8);

% Define truncation
if trunc_value ~= 0
    pts_coord = gen_truncation(pts_coord,trunc_value);
end

% Define particle dimension
x_max = max(abs(pts_coord(:,1)));
y_max = max(abs(pts_coord(:,2)));
z_max = max(abs(pts_coord(:,3)));
L_max = 2*max([x_max y_max z_max]);

% Generate mesh
sigma_particle = 0.1*L_max/10;
% sigma_particle = 0.05*ro;

% dsigma = 1E-2;
% sigma_particle = 1E-3;


sigma_o = sigma_particle;
pts_full = gen_poly(pts_coord,grid_size);
pts_full = roundn(pts_full,-3);
pts_full = unique(pts_full,'rows');
pts_full = remove_overlap(pts_full,sigma_particle,0);
while length(pts_full(:,1)) > pts_cutoff
    sigma_particle = sigma_particle + dsigma;
    pts_full = remove_overlap(pts_full,sigma_particle,0);
    disp(length(pts_full(:,1)))
end

% Minimum N
N_min = ro*sqrt(fo)*beta_kuhn;

% Determine omega
omega = sqrt(sum(pts_full.^2,2));
omega = omega/min(omega);

% Detemine shape surface area
area_scale = (sigma_o/sigma_particle)^2;
ss = alphaShape(ro*pts_coord(:,1),ro*pts_coord(:,2),ro*pts_coord(:,3));
particle_area = area_scale*surfaceArea(ss);
particle_area = 1*surfaceArea(ss);

% Find number of grafts matching grafting density
fo_array = (1:numel(pts_full(:,1)))*(0.25*pi*beta_kuhn^2)/particle_area;
fo_array = (1:numel(pts_full(:,1)))*(4*pi*beta_kuhn^2)/particle_area;
[~,indx_fo] = min(abs(fo_array-fo));
n_graft = indx_fo;

% n_graft = ceil( fo*length(pts_full(:,1))/(4*pi*beta_kuhn^2) );

disp(n_graft)
% pdf handle functions
R_chain = @(omega,N) ro*fo^(1/5)*v^(1/5)*beta_kuhn^(2/5)*N^(3/5)*omega.^(-3/5)*(beta_kuhn/ro)^(3/5) - 1.0*ro;
P_chain = @(omega,N) exp(-R_chain(omega,N).^2/(N*beta_kuhn^2) - v*fo*ro^2*N^2./((R_chain(omega,N)+ro).*omega).^3);
E_chain = @(omega,N) (-R_chain(omega,N).^2/(N*beta_kuhn^2) - v*fo*ro^2*N^2./((R_chain(omega,N)+ro).*omega).^3);
% Corona handle functions
R_chain2 = @(omega,N) ro*fo^(1/5)*v^(1/5)*beta_kuhn^(2/5)*N^(3/5)*omega.^(-3/5)*(beta_kuhn/ro)^(3/5) - 0.0*ro;
P_chain2 = @(omega,N) exp(-R_chain2(omega,N).^2/(N*beta_kuhn^2) - v*fo*ro^2*N^2./((R_chain2(omega,N)+ro).*omega).^3);

% Evaluate handle functions
P_surf = bsxfun(P_chain,omega,N);
R_surf = bsxfun(R_chain,omega,N);

% Scale probability
P_surf = (P_surf - min(P_surf(:)))/(max(P_surf(:))-min(P_surf(:)));

% Performing grafting Monte Carlo
graft_count = zeros(numel(pts_full(:,1)),1);
for ii = 1:n_mc
    fprintf('N: %i, MC Run: %i, N_min: %i\n',N,ii,int32(N_min))
    % Grafting points matrix
    graft_indx = zeros(n_graft,1);
    for j = 1:n_graft
        flag_add = 0;
        while flag_add ~= 1
            % Grafting probability
            p_jth = rand;
            % Grabbing index
            p_indx = find(P_surf >= p_jth);
            % Pick random index from selected values
            graft_jth = randi([1,length(p_indx)]);
            % Select open points
            if isempty(find(graft_indx == p_indx(graft_jth), 1))
                graft_indx(j) = p_indx(graft_jth);
                flag_add = 1;
            end
        end
    end
    % Update graft counters
    for j = 1:n_graft
        graft_count(graft_indx(j)) = graft_count(graft_indx(j)) + 1;
    end
end

% Get distances
graft_count = graft_count/n_mc;
R_graft = graft_count.*R_surf(:);

% Temp arrays
R_tmp = zeros(length(pts_full(:,1)),1);
for j = 1:numel(R_tmp(:,1))
    R_tmp(j) = (R_graft(j)+ro - min(R_graft+ro));
end
% Scale distance
R_tmp = (R_tmp - min(R_tmp(:)))/(max(R_tmp(:))-min(R_tmp(:)));
x_total = zeros(length(pts_full(:,1)),1);
y_total = zeros(length(pts_full(:,1)),1);
z_total = zeros(length(pts_full(:,1)),1);
% Generate corona
for j = 1:numel(x_total)
    n_vect = [pts_full(j,1) pts_full(j,2) pts_full(j,3)]/norm([pts_full(j,1) pts_full(j,2) pts_full(j,3)]);
    x_total(j) = n_vect(1)*R_tmp(j) + pts_full(j,1);
    y_total(j) = n_vect(2)*R_tmp(j) + pts_full(j,2);
    z_total(j) = n_vect(3)*R_tmp(j) + pts_full(j,3);
end
% Determine probablity based on new corona
omega_R = sqrt(x_total.^2+y_total.^2+z_total.^2);
omega_R = omega_R/min(omega_R);
P_surfR = bsxfun(P_chain2,omega_R,1);
P_surfR = (P_surfR - min(P_surfR(:)))/(max(P_surfR(:))-min(P_surfR(:)));

% Outputs
pts_core = pts_full;
pts_corona = [x_total y_total z_total];
pdf_core = P_surf;
pdf_corona = P_surfR;

function [fit_final,pts_poly,theta_str] = kernel_2D(n_side)

addpath('~/Documents/polyhedron_function/')

% Load in file
file_name = strcat('poly_n',num2str(n_side),'.txt');
pts_coord = textread(file_name);

% Shift angle
dangle = atan2(pts_coord(:,2),pts_coord(:,1));
dangle = sort(dangle);
indx_tmp = find(abs(dangle) > 0);
shift = min(abs(dangle(indx_tmp)));
dangle = min(diff(dangle));

% Grid size
grid_size = 150;
pts_poly = [];
% Generate grid for angular conversion
for i = 1:n_side
    % Directional vector
    if i ~= n_side
        n_v = pts_coord(i+1,:) - pts_coord(i,:);
    else
        n_v = pts_coord(1,:) - pts_coord(i,:);
    end
    lambda_space = linspace(0,norm(n_v),grid_size)';
    % Normalize
    n_v = n_v/norm(n_v);
    % Generate points
    pts_tmp = bsxfun(@times,n_v,lambda_space);
    pts_tmp = bsxfun(@plus,pts_tmp,pts_coord(i,:));
    % Store
    pts_poly = [pts_poly; pts_tmp];
end

% Reduce and rounding
pts_poly = roundn(pts_poly,-3);
pts_poly = unique(pts_poly,'rows');

% Convert to anglular kernel space
theta = atan2(pts_poly(:,2),pts_poly(:,1));

% Combine and sorting
tmp = [pts_poly theta];
tmp = sortrows(tmp,3);

% Final set of points
theta = tmp(:,3);
pts_poly = tmp(:,1:2);

% Define kernel
hr = sqrt(sum(pts_poly.^2,2));
hr = hr./min(hr);

% Angle structure
theta_str = [theta hr];

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fiting kernel function %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define inputs
vargin = cell(1,2);
vargin{1} = theta;
vargin{2} = hr;
vargin{3} = dangle;

% Fit
xo = [shift 1 1 1];
options = optimoptions('fmincon','TolX',1E-10,'Display','off');
x = fmincon(@(xo,vargin) fit_param_2D(xo,vargin),xo,[],[],[],[],[-shift 0 0 0],[shift 10 10 10],[],options,vargin);
fit_final = [x dangle];

end
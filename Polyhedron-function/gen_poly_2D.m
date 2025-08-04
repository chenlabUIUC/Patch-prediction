function pts_poly = gen_poly_2D(pts_coord,grid_size)

% Number of sides
n_side = length(pts_coord(:,1));

% Initialize
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
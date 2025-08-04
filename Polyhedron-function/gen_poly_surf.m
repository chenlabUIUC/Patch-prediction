function pts_store = gen_poly_surf(pts_core,grid_size)

% Unique points
P = unique(pts_core,'rows');

% Generate convex hull
k = boundary(P,0);

% Points generator
a1 = linspace(0,1,grid_size);
a2 = linspace(0,1,grid_size);

% Mesh grid generation
[aa1,aa2] = meshgrid(a1,a2);

% Triangulation
pts_store = cell(1,length(k(:,1)));
for kk = 1:length(k(:,1))
    % Select triangle
    triangle_test = k(kk,:);
    polygon = P(triangle_test,:);
    % Shift to origin
    polygon_shift = zeros(length(polygon(:,1)),3);
    for i = 1:length(polygon(:,1))
        polygon_shift(i,:) = polygon(i,:) - polygon(2,:);
    end
    % Leg vectors
    v1 = polygon_shift(1,:);
    v2 = polygon_shift(3,:);
    % Generate grid
    x_tmp = zeros(grid_size,grid_size);
    y_tmp = zeros(grid_size,grid_size);
    z_tmp = zeros(grid_size,grid_size);
    for i = 1:length(x_tmp(:,1))
        for j = 1:length(y_tmp(1,:))
            pts_tmp = aa1(i,j)*v1 + (1-aa1(i,j))*aa2(i,j)*v2;
            x_tmp(i,j) = pts_tmp(1) + polygon(2,1);
            y_tmp(i,j) = pts_tmp(2) + polygon(2,2);
            z_tmp(i,j) = pts_tmp(3) + polygon(2,3);
        end
    end
    % Combine arrays
    pts_final = x_tmp;
    pts_final(:,:,2) = y_tmp;
    pts_final(:,:,3) = z_tmp;
    % Store
    pts_store{kk} = pts_final;
end
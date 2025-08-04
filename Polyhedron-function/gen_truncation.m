function [pts_intersect_full] = gen_truncation(shape_coord,percent_trunc)

% Define faces
k = boundary(shape_coord,0);
% Normal vector
norm_v = sqrt(sum(shape_coord.^2,2));
n_v = bsxfun(@rdivide,shape_coord,norm_v);

% Truncation working parameters
t_range = linspace(-1,1,5000);
% Plot point and face (for truncation)
sub1 = ceil(sqrt(length(shape_coord(:,1))));
sub2 = ceil(length(shape_coord(:,1))/sub1);
pts_intersect_full = [];
for i = 1:length(shape_coord(:,1))
    % Point of truncation
    pts_trunc = n_v(i,:)*norm_v(i)*(1-percent_trunc);
    plane_trunc = [n_v(i,:) dot(n_v(i,:),-pts_trunc)];
    % Find faces with corresponding edge
    vertex_log = [];
    counter = 1;
    for j = 1:length(k(:,1))
        flag_check = sum(ismember(k(j,:),i));
        if flag_check ~= 0
            vertex_log(counter) = j;
            counter = counter + 1;
        end
    end
    indx_zr_plane = find(plane_trunc(1:3) == 0);
    pts_soln_store = cell(1,length(vertex_log));
    for j = 1:length(vertex_log)
        % Get points
        k_tmp = k(vertex_log(j),:);
        pts_tmp = shape_coord(k_tmp,:);
        % Define plane
        n1 = pts_tmp(1,:) - pts_tmp(2,:);
        n2 = pts_tmp(end,:) - pts_tmp(2,:);
        nn = cross(n2,n1);
        plane_tmp = [nn dot(nn,-pts_tmp(2,:))];
        % Define system of equations to solve
        dummy_indx = 1:3;
        indx_zr_tmp = find(plane_tmp(1:3) == 0);
        indx_zr_both = intersect(indx_zr_tmp,indx_zr_plane);
        % Finding zeros in plane normal
        if isempty(indx_zr_both)
            % If there are no common zero elements
            if (isempty(indx_zr_plane)) && (~isempty(indx_zr_tmp))
                % If truncation plane has zeros
                indx_ignore = indx_zr_tmp(1);
            elseif (isempty(indx_zr_plane)) && (isempty(indx_zr_tmp))
                indx_ignore = 1;
            else
                % If truncation plane has no zeros
                indx_ignore = indx_zr_plane(1);
            end
        else
            % If there are common zero elements
            indx_ignore = indx_zr_both(1);
        end
        % Define equations to solve
        indx_eq = find(dummy_indx ~= indx_ignore);
        mtx_tmp = [plane_trunc(indx_eq); plane_tmp(indx_eq)];
        b_tmp = -(plane_tmp(indx_ignore)*t_range + plane_tmp(4));
        b_trunc = -(plane_trunc(indx_ignore)*t_range + plane_trunc(4));
        pts_soln = zeros(length(t_range),3);
        for t = 1:length(t_range)
            soln_tmp = (mtx_tmp\[b_trunc(t);b_tmp(t)])';
            soln_store = zeros(1,3);
            indx_soln = find(dummy_indx ~= indx_ignore);
            soln_store(indx_soln) = soln_tmp;
            soln_store(indx_ignore) = t_range(t);
            pts_soln(t,:) = soln_store;
        end
        pts_soln = roundn(pts_soln,-3);
        pts_soln = unique(pts_soln,'rows');
        
        % Store for points selection
        pts_soln_store{j} = pts_soln;
    end
    
    pts_intersect = [];
    for j = 1:length(pts_soln_store)
        pts1 = pts_soln_store{j};
        for jj = j+1:length(pts_soln_store)
            pts2 = pts_soln_store{jj};
            intersect_tmp = intersect(pts1,pts2,'rows');
            if length(intersect_tmp(:,1)) < length(vertex_log)
                pts_intersect = [pts_intersect; intersect_tmp];
            end
        end
    end
    
    % Final set of points
    pts_intersect_full = [pts_intersect_full; pts_intersect];
end

pts_intersect_full = roundn(pts_intersect_full,-2);
pts_intersect_full = unique(pts_intersect_full,'rows');
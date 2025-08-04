function [pts_store,pts_pre,pts_indv_store] = vertex_truncation(pts_core,frac_trunc)

% Define vertex connectivity
vertex_shortest = zeros(1,length(pts_core(:,1)));
vertex_connect = cell(1,length(pts_core(:,1)));
for i = 1:length(pts_core(:,1))
    % Distance
    dr_tmp = bsxfun(@minus,pts_core,pts_core(i,:));
    dr_tmp = sqrt(sum(dr_tmp.^2,2));    
    % Remove zero point
    indx_sel = find(dr_tmp ~= 0);
    dr_sel = dr_tmp(indx_sel);    
    % Sort
    dr_sort = sort(dr_sel);
    % Select closest points
    indx_sel2 = find((dr_sel-dr_sort(3)) < 1E-2);
    indx_sel_final = indx_sel(indx_sel2);
    % Shortest distance
    vertex_shortest(i) = min(dr_sel);
    % Store
    vertex_connect{i} = indx_sel_final;
end

% Truncate points
pts_store = [];
pts_pre = [];
pts_indv_store = cell(1,length(pts_core(:,1)));
for i = 1:length(pts_core(:,1))
    % Get point
    pts_tmp = pts_core(i,:);
    % Get shortened point
    nv_tmp = pts_tmp/norm(pts_tmp);
    pts_tmp_t = pts_tmp - frac_trunc*norm(pts_tmp)*nv_tmp;
    pts_pre = [pts_pre; pts_tmp_t];
    % Get project vertex
    indx_project = vertex_connect{i};
    % Project points onto line
    pts_indv = [];
    for j = 1:length(indx_project)       
        % Determine distance        
        lo = pts_tmp;
        lp = (pts_core(indx_project(j),:)-pts_tmp)/norm(pts_core(indx_project(j),:)-pts_tmp);
        n = (pts_tmp-pts_tmp_t)/norm(pts_tmp-pts_tmp_t);
        pt = pts_tmp_t;
        dmove = (dot(n,pt)-dot(n,lo))/(dot(n,lp));
        % Determine how far to move
        prj_tmp = pts_tmp + lp*dmove;
        pts_store = [pts_store; prj_tmp];
        % Store subset
        pts_indv = [pts_indv; prj_tmp];
    end
    pts_indv_store{i} = pts_indv;
end

end
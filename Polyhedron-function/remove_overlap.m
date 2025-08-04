function [pts_overlap] = remove_overlap(pts_poly,sigma,percent_overlap)
% Randomize data
for i = 1:5E3
    pts_poly = pts_poly(randperm(size(pts_poly,1)),:);
end
% Neighbor list
neigh_cell = cell(1,length(pts_poly(:,1)));
% fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
% fprintf('Removing Overlaping Particles\n')
% fprintf('Building Cell List...')
dr2 = pdist2(pts_poly,pts_poly);
for i = 1:length(pts_poly(:,1))
    % Get distance
%     dr = bsxfun(@minus,pts_poly,pts_poly(i,:));
%     norm_dr = sqrt(sum(dr.^2,2));
    norm_dr = dr2(i,:);
    % Remove central particle
    indx_cell = find(norm_dr < (1-percent_overlap)*sigma);
    indx_cell = setdiff(indx_cell,i);
    % Store neighbor list
    neigh_cell{i} = indx_cell;
end
% fprintf('done\n')
% fprintf('Remove particle...')
% Removal check
pts_overlap = pts_poly(1,:);
indx_particle = 1;
% Loop through particles
for i = 2:1:length(pts_poly(:,1))
    % Loop through kept points
    indx_ref = neigh_cell{i};
    indx_check = intersect(indx_ref,indx_particle);
    indx_flag = 0;
    if ~isempty(indx_check)
        indx_flag = 1;
    end
    if indx_flag == 0
        indx_particle(end+1) = i;
        pts_overlap = [pts_overlap; pts_poly(i,:)];
    end
end
% fprintf('done\n')
% fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
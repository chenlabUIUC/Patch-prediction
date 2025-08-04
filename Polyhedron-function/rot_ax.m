function rtn_mtx = rot_ax(dr,flag_z)
new_z = flag_z*dr;
new_z = new_z/norm(new_z);
new_y = cross([1 0 0],new_z);
if norm(new_y) == 0
    new_y = [0 sign(new_z(1)) 0];
end
new_x = cross(new_z,new_y);

new_x = new_x/norm(new_x);
new_y = new_y/norm(new_y);
new_z = new_z/norm(new_z);
% Put into matrix
rtn_mtx = [new_z; new_y; new_x]';
end

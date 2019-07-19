function [cells, vertices] = convert_sgrid_to_hex_ugrid(is_cell_active,is_vertex_active,xv,yv,zv)


[nvx,nvy,nvz] = size(xv);

vertex_ids = compute_active_vertex_ids(is_vertex_active);

% determine total number of active hexs
ncells = sum(sum(sum(is_cell_active)));

% determine the subscript of active hexs
[i1,i2,i3]=ind2sub(size(is_cell_active),find(is_cell_active==1));


% set up the vertex ids for all active hex
cells = [ 8*ones(ncells,1) ...
    vertex_ids(sub2ind([nvx nvy nvz],i1  ,i2  ,i3  )) ...
    vertex_ids(sub2ind([nvx nvy nvz],i1+1,i2  ,i3  )) ...
    vertex_ids(sub2ind([nvx nvy nvz],i1+1,i2+1,i3  )) ...
    vertex_ids(sub2ind([nvx nvy nvz],i1  ,i2+1,i3  )) ...
    vertex_ids(sub2ind([nvx nvy nvz],i1  ,i2  ,i3+1)) ...
    vertex_ids(sub2ind([nvx nvy nvz],i1+1,i2  ,i3+1)) ...
    vertex_ids(sub2ind([nvx nvy nvz],i1+1,i2+1,i3+1)) ...
    vertex_ids(sub2ind([nvx nvy nvz],i1  ,i2+1,i3+1)) ];

% set the (x,y,z) of active vertices
loc = find(is_vertex_active == 1);

vertices = [xv(loc) yv(loc) zv(loc)];

end

function is_vertex_active = identify_active_vertices(sgrid, is_cell_active)

nx = sgrid.nx;
ny = sgrid.ny;
nz = sgrid.nz;

nvx = nx + 1;
nvy = ny + 1;
nvz = nz + 1;

is_vertex_active = zeros(nvx,nvy,nvz);

[i1,i2,i3]=ind2sub(size(is_cell_active),find(is_cell_active==1));

is_vertex_active(sub2ind([nvx nvy nvz],i1  ,i2  ,i3  )) = 1;
is_vertex_active(sub2ind([nvx nvy nvz],i1+1,i2  ,i3  )) = 1;
is_vertex_active(sub2ind([nvx nvy nvz],i1+1,i2+1,i3  )) = 1;
is_vertex_active(sub2ind([nvx nvy nvz],i1  ,i2+1,i3  )) = 1;
is_vertex_active(sub2ind([nvx nvy nvz],i1  ,i2  ,i3+1)) = 1;
is_vertex_active(sub2ind([nvx nvy nvz],i1+1,i2  ,i3+1)) = 1;
is_vertex_active(sub2ind([nvx nvy nvz],i1+1,i2+1,i3+1)) = 1;
is_vertex_active(sub2ind([nvx nvy nvz],i1  ,i2+1,i3+1)) = 1;


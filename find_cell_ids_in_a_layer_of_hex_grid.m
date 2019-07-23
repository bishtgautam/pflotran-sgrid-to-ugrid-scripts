function cell_idx = find_cell_ids_in_a_layer_of_hex_grid(sgrid,cell_ids,layer_offset_from_surface)

nx = sgrid.nx;
ny = sgrid.ny;

cell_idx = zeros(nx*ny,1);
count = 0;
for jj = 1:ny
    for ii = 1:nx
        count = count + 1;
        loc = find(cell_ids(ii,jj,:)>0);
        cell_idx(count) = cell_ids(ii,jj,loc(end-layer_offset_from_surface));
    end
end

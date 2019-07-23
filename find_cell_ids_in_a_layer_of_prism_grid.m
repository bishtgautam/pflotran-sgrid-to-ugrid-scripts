function cell_idx = find_cell_ids_in_a_layer_of_prism_grid(sgrid,nz_prism,layer_offset_from_surface)

nx = sgrid.nx;
ny = sgrid.ny;


idx_beg = 4*nx*ny*(nz_prism-layer_offset_from_surface-1)+1;
idx_end = idx_beg + 4*nx*ny -1;

cell_idx = [idx_beg:idx_end]';



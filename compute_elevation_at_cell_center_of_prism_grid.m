function zc_top = compute_elevation_at_cell_center_of_prism_grid(sgrid,is_cell_active,z)

nx = sgrid.nx;
ny = sgrid.ny;

zc_top = reshape(z(reshape(sum(is_cell_active,3),nx*ny,1)),nx,ny);

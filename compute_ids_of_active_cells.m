function cell_ids = compute_ids_of_active_cells(sgrid, h5_material_filename)


nx = sgrid.nx;
ny = sgrid.ny;
nz = sgrid.nz;

cell_ids   = zeros(nx,ny,nz);

if (~isempty(h5_material_filename))
    mat_ids = h5read(h5_material_filename,'/Materials/Material Ids');
    loc = find(mat_ids >  0);
    
    cell_ids(loc) = [1:length(loc)];
else
    loc = [1:nx*ny*nz];
    cell_ids(loc) = loc;
end

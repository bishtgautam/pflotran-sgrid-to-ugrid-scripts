function is_cell_active = identify_active_cells(sgrid, h5_material_filename)

nx = sgrid.nx;
ny = sgrid.ny;
nz = sgrid.nz;

is_cell_active   = ones(nx,ny,nz);

if (~isempty(h5_material_filename))
    mat_ids = h5read(h5_material_filename,'/Materials/Material Ids');
    
    loc = find(mat_ids == 0);
    
    is_cell_active(loc) = 0;
end

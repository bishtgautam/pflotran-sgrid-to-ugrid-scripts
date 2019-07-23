function [ugrid_mat_ids, ugrid_mat_cell_ids] = compute_ugrid_materials(sgrid, h5_material_filename, cell_ids)

nx = sgrid.nx;
ny = sgrid.ny;
nz = sgrid.nz;


if (~isempty(h5_material_filename))
    mat_ids      = h5read(h5_material_filename,'/Materials/Material Ids');
    mat_cell_ids = h5read(h5_material_filename,'/Materials/Cell Ids');
    
    if (nx*ny*nz ~= length(mat_ids))
        error(sprintf("The number of grid cells are not equal to material ids.\nNo. of grid cells   = %d\nNo. of material ids = %d%d",nx*ny*nz,length(mat_ids)))
    end
    
    loc = find(mat_ids >  0);
    ugrid_mat_ids = mat_ids(loc);
    ugrid_mat_cell_ids = cell_ids(mat_cell_ids(loc));
else
    ugrid_mat_ids      = ones(nx*ny*nz,1);
    ugrid_mat_cell_ids = [1:nx*ny*nz];
end
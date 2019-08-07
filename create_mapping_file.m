function create_mapping_file(sgrid,ugrid,elm_nsoil)

MESH_PRI = 1;
MESH_HEX = 2;


h5_ugrid_filename = ugrid.h5_filename;
mesh_typename     = ugrid.mesh_typename;
map_base_filename = ugrid.map_base_filename;

switch lower(mesh_typename)
    case 'prism'
        mesh_type = MESH_PRI;
    case 'hex'
        mesh_type = MESH_HEX;
    otherwise
        error(['Unknown mesh_typename: ' mesh_typename]);
end

nx = sgrid.nx;
ny = sgrid.ny;
nz = sgrid.nz;

top = h5read(h5_ugrid_filename,'/Regions/Top/Vertex Ids')';
ncells_per_lyr = size(top,1);

h5_material_filename = sgrid.h5_material_filename;

cell_ids = compute_ids_of_active_cells(sgrid, h5_material_filename);


switch mesh_type
    case MESH_HEX
        pfl_cell_idx = zeros(ncells_per_lyr,elm_nsoil);
        elm_cell_idx = zeros(ncells_per_lyr,elm_nsoil);
        pfl_2_elm_wt = ones(ncells_per_lyr,elm_nsoil);
        elm_2_pfl_wt = ones(ncells_per_lyr,elm_nsoil);

        for kk = 1:elm_nsoil
            pfl_cell_idx(:,kk) = find_cell_ids_in_a_layer_of_hex_grid(sgrid,cell_ids,kk-1);
            elm_cell_idx(:,kk) = [kk:elm_nsoil:nx*ny*elm_nsoil];
        end
        
    case MESH_PRI
        if (~isfield(ugrid,'nz'))
            zc_top   = compute_elevation_at_cell_center_of_prism_grid(sgrid,is_cell_active,z);
            zv_top   = compute_elevation_at_vertices_of_prism_grid(sgrid, zc_top);
            nz_prism = ceil((max(max(zv_top))-z_min)/dz);
        else
            nz_prism = ugrid.nz;
        end
        
        pfl_cell_idx = zeros(ncells_per_lyr,elm_nsoil);
        elm_cell_idx = zeros(ncells_per_lyr,elm_nsoil);
        pfl_2_elm_wt = ones(ncells_per_lyr,elm_nsoil)/4.0;
        elm_2_pfl_wt = ones(ncells_per_lyr,elm_nsoil);

        for kk = 1:elm_nsoil
            pfl_cell_idx(:,kk) = find_cell_ids_in_a_layer_of_prism_grid(sgrid,nz_prism,kk-1);
            elm_cell_idx(1:4:end,kk) = [kk:elm_nsoil:nx*ny*elm_nsoil];
            elm_cell_idx(2:4:end,kk) = [kk:elm_nsoil:nx*ny*elm_nsoil];
            elm_cell_idx(3:4:end,kk) = [kk:elm_nsoil:nx*ny*elm_nsoil];
            elm_cell_idx(4:4:end,kk) = [kk:elm_nsoil:nx*ny*elm_nsoil];
        end
end

[a,b]=size(elm_cell_idx);
elm_cell_idx_1d = reshape(elm_cell_idx,a*b,1);
pfl_cell_idx_1d = reshape(pfl_cell_idx,a*b,1);
pfl_2_elm_wt_1d = reshape(pfl_2_elm_wt,a*b,1);
elm_2_pfl_wt_1d = reshape(elm_2_pfl_wt,a*b,1);


% PF -- to -- ELM subsurface
[~,idx] = sort(elm_cell_idx_1d);
map(:,1) = elm_cell_idx_1d(idx);
map(:,2) = pfl_cell_idx_1d(idx);
map(:,3) = pfl_2_elm_wt_1d(idx);
write_mapping_file([map_base_filename '_pf2elm.meshmap'], elm_nsoil, map);
clear map

% ELM --to -- PF subsurface
[~,idx] = sort(pfl_cell_idx_1d);
map(:,1) = pfl_cell_idx_1d(idx);
map(:,2) = elm_cell_idx_1d(idx);
map(:,3) = elm_2_pfl_wt_1d(idx);
write_mapping_file([map_base_filename '_elm2pf.meshmap'], elm_nsoil, map);
clear map

% ELM --to -- PF surface only
[~,idx] = sort(pfl_cell_idx(:,1));
map(:,1) = pfl_cell_idx(idx,1);
map(:,2) = elm_cell_idx(idx,1);
map(:,3) = elm_2_pfl_wt(idx,1);
write_mapping_file([map_base_filename '_elm2pf_surf.meshmap'], elm_nsoil, map);
clear map

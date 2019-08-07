function [ugrid_mat_ids, ugrid_mat_cell_ids] = compute_prism_ugrid_materials(sgrid, h5_material_filename, zv_prism)

nx = sgrid.nx;
ny = sgrid.ny;
nz = sgrid.nz;

dz = sgrid.dz;

nz_prism = length(zv_prism) - 1;

if (~isempty(h5_material_filename))
    mat_ids      = h5read(h5_material_filename,'/Materials/Material Ids');
    mat_cell_ids = h5read(h5_material_filename,'/Materials/Cell Ids');
    
    mids(mat_cell_ids) = mat_ids;
    
    clear mat_ids
    
    mids = reshape(mids,nx,ny,nz);
    
    for jj = 1:ny
        for ii = 1:nx
            
            % Going up from the bottom, determine the index to the last
            % active layer
            kk_max = max(find(mids(ii,jj,:)>0));
            
            % Going up from the bottom, extract the material ids
            tmp_ids = reshape(mids(ii,jj,1:kk_max),kk_max,1);
            
            % Now flip the material ids, so we will be traversing from the top
            % layer down
            tmp_ids = flipud(tmp_ids);
            
            zc = [dz/2:dz:dz*kk_max];
            
            for kk = 1:nz_prism
                loc = find(zc>= zv_prism(kk) & zc<=zv_prism(kk+1));
                
                if (~isempty(loc))
                    mat_ids(ii,jj,nz_prism-kk+1) = mode(tmp_ids(loc));
                else
                    mat_ids(ii,jj,nz_prism-kk+1) = tmp_ids(end);
                end
                
            end
        end
    end
    
    tmp = reshape(mat_ids,nx*ny*nz_prism,1);
    ugrid_mat_ids = reshape([tmp tmp tmp tmp]',nx*ny*nz_prism*4,1);
    ugrid_mat_cell_ids = [1:nx*ny*nz_prism*4];
    
    
else
    ugrid_mat_ids      = ones(nx*ny*nz_prism*4,1);
    ugrid_mat_cell_ids = [1:nx*ny*nz_prism*5];
end
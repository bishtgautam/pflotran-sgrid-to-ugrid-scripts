function convert_sgrid_to_prism_ugrid_driver(sgrid,h5_ugrid_filename)

nx = sgrid.nx;
ny = sgrid.ny;
nz = sgrid.nz;
dx = sgrid.dx;
dy = sgrid.dy;
dz = sgrid.dz;

h5_material_filename = sgrid.h5_material_filename;
h5_region_filename   = sgrid.h5_region_filename;


x_min = sgrid.origin_x;
y_min = sgrid.origin_y;
z_min = sgrid.origin_z;


mat_ids = h5read(h5_material_filename,'/Materials/Material Ids');
mat_cell_ids = h5read(h5_material_filename,'/Materials/Cell Ids');

if (nx*ny*nz ~= length(mat_ids))
    error(sprintf("The number of grid cells are not equal to material ids.\nNo. of grid cells   = %d\nNo. of material ids = %d%d",nx*ny*nz,length(mat_ids)))
end


x = [0:dx:dx*nx] + x_min;
y = [0:dx:dy*ny] + y_min;
z = [0:dz:dz*nz] + z_min;

nvx = nx + 1;
nvy = ny + 1;
nvz = nz + 1;

xv               = zeros(nvx,nvy,nvz);
yv               = zeros(nvx,nvy,nvz);
zv               = zeros(nvx,nvy,nvz);

is_cell_active   = identify_active_cells(sgrid, h5_material_filename);
cell_ids         = compute_ids_of_active_cells(sgrid, h5_material_filename);
%is_vertex_active = identify_active_vertices(sgrid, is_cell_active);

%[ugrid_mat_ids, ugrid_mat_cell_ids] = compute_ugrid_materials(sgrid, h5_material_filename, cell_ids);

for kk = 1:nvz
    for jj = 1:nvy
        for ii = 1:nvx
            xv(ii,jj,kk) = x(ii);
            yv(ii,jj,kk) = y(jj);
            zv(ii,jj,kk) = z(kk);
        end
    end
end

zc_top   = compute_elevation_at_cell_center_of_prism_grid(sgrid,is_cell_active,z);
zv_top   = compute_elevation_at_vertices_of_prism_grid(sgrid, zc_top);

if (~isfield(sgrid,'nz_prism'))
    nz_prism = ceil((max(max(zv_top))-z_min)/dz);
else
    nz_prism = sgrid.nz_prism;
end

[cells, vertices] = convert_sgrid_to_prism_ugrid(xv(:,:,1),yv(:,:,1),zv_top,nz_prism,dz);

ugrid_mat_ids = ones(size(cells,1),1);
ugrid_mat_cell_ids = [1:size(cells,1)];


system(['rm -f ' h5_ugrid_filename]);

fprintf('Creating the following mesh file:\n\t%s\n\n', h5_ugrid_filename);

h5create(h5_ugrid_filename,'/Domain/Cells',size(cells'),'Datatype','int64');
h5create(h5_ugrid_filename,'/Domain/Vertices',size(vertices'));
h5create(h5_ugrid_filename,'/Materials/Cell Ids',length(ugrid_mat_ids),'Datatype','int64');
h5create(h5_ugrid_filename,'/Materials/Material Ids',length(ugrid_mat_cell_ids),'Datatype','int64');


h5write(h5_ugrid_filename,'/Domain/Cells',int64(cells'));
h5write(h5_ugrid_filename,'/Domain/Vertices',vertices');
h5write(h5_ugrid_filename,'/Materials/Cell Ids',int64(ugrid_mat_cell_ids));
h5write(h5_ugrid_filename,'/Materials/Material Ids',int64(ugrid_mat_ids));

top_cell_idx = find_cell_ids_in_a_layer_of_prism_grid(sgrid,nz_prism,0);
top_cells = cells(top_cell_idx,:);

region_info = h5info(h5_region_filename,'/Regions');

for rr = 1:length(region_info.Groups)

    if (rr == 1)
        disp('Adding following regions: ');
    end
    disp([' ' region_info.Groups(rr).Name]);
    info = h5info(h5_region_filename,region_info.Groups(rr).Name);
    cids = h5read(h5_region_filename,[region_info.Groups(rr).Name '/' info.Datasets(1).Name]);
    
    if (~strcmp(info.Datasets(2).Name,'Face Ids'))
        error(['For ' region_info.Groups(rr).Name ': "Face Ids" not found']);
    end
    
    fids = h5read(h5_region_filename,[region_info.Groups(rr).Name '/' info.Datasets(2).Name]);
    

    cids = cids(find(fids == 6));
    ugrid_region_fids = zeros(length(cids)*4,3);

    tmp_cids = double(cids) - floor(double(cids)/nx/ny)*nx*ny;

    h5create(h5_ugrid_filename,[region_info.Groups(rr).Name '/Vertex Ids'],size(ugrid_region_fids'),'Datatype','int64');

    for ii = 1:length(cids)
        
        tmp_id = cids(ii) - floor(double(cids(ii))/nx/ny)*nx*ny;
        
        ugrid_region_fids((ii-1)*4+1:ii*4,:) = top_cells((tmp_id-1)*4+1:tmp_id*4,5:7);
    end
    
    h5write(h5_ugrid_filename,[region_info.Groups(rr).Name '/Vertex Ids'],int64(ugrid_region_fids'));

end

disp(' /Regions/All')
h5create(h5_ugrid_filename,['/Regions/All/Cell Ids'],length(ugrid_mat_cell_ids),'Datatype','int64');
h5write(h5_ugrid_filename,['/Regions/All/Cell Ids'],int64(ugrid_mat_cell_ids));

disp(' /Regions/Top')
h5create(h5_ugrid_filename,['/Regions/Top/Vertex Ids'],size(top_cells(:,5:7)'),'Datatype','int64');
h5write(h5_ugrid_filename,['/Regions/Top/Vertex Ids'],int64(top_cells(:,5:7)'));

end
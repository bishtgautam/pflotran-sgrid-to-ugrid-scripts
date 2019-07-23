function convert_sgrid_to_ugrid(sgrid,h5_ugrid_filename)

nx = sgrid.nx;
ny = sgrid.ny;
nz = sgrid.nz;
dx = sgrid.dx;
dy = sgrid.dy;
dz = sgrid.dz;

h5_material_filename = '';
h5_region_filename   = '';

if isfield(sgrid,'h5_material_filename')
 h5_material_filename = sgrid.h5_material_filename;
end

if isfield(sgrid,'h5_region_filename')
    h5_region_filename   = sgrid.h5_region_filename;
end

x_min = sgrid.origin_x;
y_min = sgrid.origin_y;
z_min = sgrid.origin_z;

x = [0:dx:dx*nx] + x_min;
y = [0:dx:dy*ny] + y_min;
z = [0:dz:dz*nz] + z_min;

nvx = nx + 1;
nvy = ny + 1;
nvz = nz + 1;

is_cell_active   = identify_active_cells(sgrid, h5_material_filename);
cell_ids         = compute_ids_of_active_cells(sgrid, h5_material_filename);
is_vertex_active = identify_active_vertices(sgrid, is_cell_active);

[ugrid_mat_ids, ugrid_mat_cell_ids] = compute_ugrid_materials(sgrid, h5_material_filename, cell_ids);

xv = zeros(nvx,nvy,nvz);
yv = zeros(nvx,nvy,nvz);
zv = zeros(nvx,nvy,nvz);

for kk = 1:nvz
    for jj = 1:nvy
        for ii = 1:nvx
            xv(ii,jj,kk) = x(ii);
            yv(ii,jj,kk) = y(jj);
            zv(ii,jj,kk) = z(kk);
        end
    end
end

[cells, vertices] = convert_sgrid_to_hex_ugrid(is_cell_active,is_vertex_active,xv,yv,zv);

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


if (~isempty(h5_region_filename))
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
        
        
        ugrid_region_cids   = cell_ids(cids);
        tmp_cell_vertex_ids = cells(cell_ids(cids),2:end);
        ugrid_region_fids = zeros(length(cids),4);
        
        h5create(h5_ugrid_filename,[region_info.Groups(rr).Name '/Vertex Ids'],size(ugrid_region_fids'),'Datatype','int64');
        
        for ii = 1:length(cids)
            switch fids(ii)
                case 1 % west
                    tmp_face_ids = [1 5 8 4];
                case 2 % east
                    tmp_face_ids = [2 3 7 6];
                case 3 % south
                    tmp_face_ids = [1 2 6 5];
                case 4 % north
                    tmp_face_ids = [3 4 8 7];
                case 5 % bottom
                    tmp_face_ids = [1 4 3 2];
                case 6 % top
                    tmp_face_ids = [5 6 7 8];
                otherwise
                    error('Invalid face id')
            end
            ugrid_region_fids(ii,:) = tmp_cell_vertex_ids(ii,tmp_face_ids);
        end
        
        h5write(h5_ugrid_filename,[region_info.Groups(rr).Name '/Vertex Ids'],int64(ugrid_region_fids'));
        
    end
end

cids = reshape(sum(is_cell_active,3),nx*ny,1);

top_cell_idx = find_cell_ids_in_a_layer_of_hex_grid(sgrid,cell_ids,0);

disp(' /Regions/Top');
ugrid_region_fids = cells(top_cell_idx,6:end);
h5create(h5_ugrid_filename,'/Regions/Top/Vertex Ids',size(ugrid_region_fids'),'Datatype','int64');
h5write(h5_ugrid_filename,'/Regions/Top/Vertex Ids',int64(ugrid_region_fids'));

end
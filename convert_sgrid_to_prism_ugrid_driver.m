function convert_sgrid_to_prism_ugrid_driver(sgrid,h5_material_filename,h5_region_filename,h5_ugrid_filename)

nx = sgrid.nx;
ny = sgrid.ny;
nz = sgrid.nz;
dx = sgrid.dx;
dy = sgrid.dy;
dz = sgrid.dz;
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

is_cell_active   = ones(nx,ny,nz);
cell_ids         = zeros(nx,ny,nz);
xv               = zeros(nvx,nvy,nvz);
yv               = zeros(nvx,nvy,nvz);
zv               = zeros(nvx,nvy,nvz);
is_vertex_active = zeros(nvx,nvy,nvz);


loc = find(mat_ids == 0);is_cell_active(loc) = 0;
loc = find(mat_ids >  0);cell_ids(loc) = [1:length(loc)];
ugrid_mat_ids = mat_ids(loc);
ugrid_mat_cell_ids = cell_ids(mat_cell_ids(loc));

[i1,i2,i3]=ind2sub(size(is_cell_active),find(is_cell_active==1));

is_vertex_active(sub2ind([nvx nvy nvz],i1  ,i2  ,i3  )) = 1;
is_vertex_active(sub2ind([nvx nvy nvz],i1+1,i2  ,i3  )) = 1;
is_vertex_active(sub2ind([nvx nvy nvz],i1+1,i2+1,i3  )) = 1;
is_vertex_active(sub2ind([nvx nvy nvz],i1  ,i2+1,i3  )) = 1;
is_vertex_active(sub2ind([nvx nvy nvz],i1  ,i2  ,i3+1)) = 1;
is_vertex_active(sub2ind([nvx nvy nvz],i1+1,i2  ,i3+1)) = 1;
is_vertex_active(sub2ind([nvx nvy nvz],i1+1,i2+1,i3+1)) = 1;
is_vertex_active(sub2ind([nvx nvy nvz],i1  ,i2+1,i3+1)) = 1;


for kk = 1:nvz
    for jj = 1:nvy
        for ii = 1:nvx
            xv(ii,jj,kk) = x(ii);
            yv(ii,jj,kk) = y(jj);
            zv(ii,jj,kk) = z(kk);
        end
    end
end

zc_top = reshape(z(reshape(sum(is_cell_active,3),nx*ny,1)),nx,ny);

zv_top = zv(:,:,1)*0;

zv_top(2:nx,2:ny) = ...
    (...
    zc_top(1:nx-1,1:ny-1) + ...
    zc_top(2:nx  ,1:ny-1) + ...
    zc_top(1:nx-1,2:ny  ) + ...
    zc_top(2:nx  ,2:ny  ) ...
    )/4;

zv_top(2:nx,1   ) = (zc_top(1:nx-1,1     ) + zc_top(2:nx  ,1   ))/2;
zv_top(2:nx,ny+1  ) = (zc_top(1:nx-1,ny    ) + zc_top(2:nx  ,ny  ))/2;
zv_top(1   ,2:ny) = (zc_top(1     ,1:ny-1) + zc_top(1     ,2:ny))/2;
zv_top(nx+1  ,2:ny) = (zc_top(nx    ,1:ny-1) + zc_top(nx    ,2:ny))/2;

zv_top(1,1) = zc_top(1,1);
zv_top(1,ny+1) = zc_top(1,ny);
zv_top(nx+1,1) = zc_top(nx,1);
zv_top(nx+1,ny+1) = zc_top(nx,ny);


%zv_top = reshape(z(reshape(sum(is_vertex_active,3),nvx*nvy,1)),nvx,nvy);

[cells, vertices] = convert_sgrid_to_prism_ugrid(xv(:,:,1),yv(:,:,1),zv_top,ceil((max(max(zv_top))-z_min)/dz),dz);

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

top_cells = cells(end-nx*ny*4+1:end,:);



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
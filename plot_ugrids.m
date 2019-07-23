function plot_ugrids(save_figure)

h5_ugrid_filename = '/Users/bish218/projects_1drive/elm-pflotran-coupling/dataset/5km/pflotran/5km_hex_ugrid.h5';
region_info = h5info(h5_ugrid_filename,'/Regions');

vertices = h5read(h5_ugrid_filename,'/Domain/Vertices')';

figure;
for rr = 1:length(region_info.Groups)
    if (~strcmp(region_info.Groups(rr).Name,'/Regions/Top'))
        info = h5info(h5_ugrid_filename,region_info.Groups(rr).Name);
        cids = h5read(h5_ugrid_filename,[region_info.Groups(rr).Name '/' info.Datasets(1).Name]);

        if (~strcmp(info.Datasets(1).Name,'Vertex Ids'))
            error(['For ' region_info.Groups(rr).Name ': "Vertex Ids" not found']);
        end

        vids = h5read(h5_ugrid_filename,[region_info.Groups(rr).Name '/' info.Datasets(1).Name])';

        for k = 1:size(vids,1)
            xx = vertices([vids(k,:) vids(k,1) ],1); yy = vertices([vids(k,:) vids(k,1) ],2); zz = vertices([vids(k,:) vids(k,1) ],3);
            plot3(xx,yy,zz,'-rs','linewidth',2); hold all;
        end
    end
end

vids = h5read(h5_ugrid_filename,'/Regions/Top/Vertex Ids')';
for k = 1:size(vids,1)
    xx = vertices([vids(k,:) vids(k,1) ],1); yy = vertices([vids(k,:) vids(k,1) ],2); zz = vertices([vids(k,:) vids(k,1) ],3);
    %plot3(xx,yy,zz,'-bs','linewidth',2); hold all;
    fill3(xx,yy,zz,mean(zz)); hold all;
end
plot3(vertices(:,1),vertices(:,2),vertices(:,3),'.k')
grid on
zlim([0 800])
view(-33,70)
caxis([100 500])
colorbar
if (save_figure)
    orient landscape
    print -dpdf hex_ugrid.pdf -fillpage
end


h5_ugrid_filename = '/Users/bish218/projects_1drive/elm-pflotran-coupling/dataset/5km/pflotran/5km_prism_ugrid.h5';
vertices = h5read(h5_ugrid_filename,'/Domain/Vertices')';
cells = h5read(h5_ugrid_filename,'/Domain/Cells')';

figure;
vids = h5read(h5_ugrid_filename,'/Regions/Top/Vertex Ids')';
for k = 1:size(vids,1)
    xx = vertices([vids(k,:) vids(k,1) ],1); yy = vertices([vids(k,:) vids(k,1) ],2); zz = vertices([vids(k,:) vids(k,1) ],3);
    plot3(xx,yy,zz,'-ko','linewidth',2); hold all;
end

% figure;
vids = h5read(h5_ugrid_filename,'/Regions/Top/Vertex Ids')';
for k = 1:size(vids,1)
    xx = vertices([vids(k,:) vids(k,1) ],1); yy = vertices([vids(k,:) vids(k,1) ],2); zz = vertices([vids(k,:) vids(k,1) ],3);
    fill3(xx(1:3),yy(1:3),zz(1:3),mean(zz(1:3))); hold all;
end

region_info = h5info(h5_ugrid_filename,'/Regions');
for rr = 1:length(region_info.Groups)
    if (~strcmp(region_info.Groups(rr).Name,'/Regions/Top'))
    info = h5info(h5_ugrid_filename,region_info.Groups(rr).Name);
    cids = h5read(h5_ugrid_filename,[region_info.Groups(rr).Name '/' info.Datasets(1).Name]);
    
    if (strcmp(info.Datasets(1).Name,'Vertex Ids'))
        
        vids = h5read(h5_ugrid_filename,[region_info.Groups(rr).Name '/' info.Datasets(1).Name])';
        
        for k = 1:size(vids,1)
            xx = vertices([vids(k,:) vids(k,1) ],1); yy = vertices([vids(k,:) vids(k,1) ],2); zz = vertices([vids(k,:) vids(k,1) ],3);
            plot3(xx,yy,zz,'-rs','linewidth',2); hold all;
        end
    end
    end
end
grid on;
zlim([0 800])
view(-33,70)
caxis([100 500])
colorbar
if (save_figure)
    orient landscape
    print -dpdf prism_ugrid.pdf -fillpage
end


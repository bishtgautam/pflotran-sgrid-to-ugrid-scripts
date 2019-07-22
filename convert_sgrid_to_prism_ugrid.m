function [cells, vertices] = convert_sgrid_to_prism_ugrid(xv,yv,zv,nz,dz)

[nvx,nvy] = size(xv);

nx = nvx-1;
ny = nvy-1;

nvz = nz + 1;

vertices = zeros(nvx*nvy*nvz + nx*ny*nvz,3);
cells    = zeros(nx *ny *nz * 4, 7);

id_cen = zeros(nx ,ny ,nvz);



id_v = reshape([1:nvx*nvy*nvz],nvx,nvy,nvz);

count = 1;
for kk = 1:nvz
    for jj = 1:nvy
        for ii = 1:nvx
            vertices(count,1) = xv(ii,jj);
            vertices(count,2) = yv(ii,jj);
            vertices(count,3) = zv(ii,jj) - (nvz-kk)*dz;
            count = count + 1;
        end
    end
end

for kk = 1:nvz
    for jj = 1:ny
        for ii = 1:nx
            vertices(count,1) = (...
                xv(ii  ,jj  ) + ...
                xv(ii+1,jj  ) + ...
                xv(ii+1,jj+1) + ...
                xv(ii  ,jj+1) )/4;
            
            vertices(count,2) = (...
                yv(ii  ,jj  ) + ...
                yv(ii+1,jj  ) + ...
                yv(ii+1,jj+1) + ...
                yv(ii  ,jj+1) )/4;
            
            vertices(count,3) = (...
                zv(ii  ,jj  ) + ...
                zv(ii+1,jj  ) + ...
                zv(ii+1,jj+1) + ...
                zv(ii  ,jj+1) )/4 - (nvz-kk)*dz;
            
            id_cen(ii,jj,kk) = count;
            count = count + 1;
        end
    end
end

cells(:,1) = 6;
count = 0;
for kk = 1:nz
    for jj = 1:ny
        for ii = 1:nx
            count = count + 1; cells(count,2:7) = [id_v(ii  ,jj  ,kk  ) id_v(ii+1,jj  ,kk  ) id_cen(ii,jj,kk)  id_v(ii  ,jj  ,kk+1) id_v(ii+1,jj  ,kk+1) id_cen(ii,jj,kk+1)];
            count = count + 1; cells(count,2:7) = [id_v(ii+1,jj  ,kk  ) id_v(ii+1,jj+1,kk  ) id_cen(ii,jj,kk)  id_v(ii+1,jj  ,kk+1) id_v(ii+1,jj+1,kk+1) id_cen(ii,jj,kk+1)];
            count = count + 1; cells(count,2:7) = [id_v(ii+1,jj+1,kk  ) id_v(ii  ,jj+1,kk  ) id_cen(ii,jj,kk)  id_v(ii+1,jj+1,kk+1) id_v(ii  ,jj+1,kk+1) id_cen(ii,jj,kk+1)];
            count = count + 1; cells(count,2:7) = [id_v(ii  ,jj+1,kk  ) id_v(ii  ,jj  ,kk  ) id_cen(ii,jj,kk)  id_v(ii  ,jj+1,kk+1) id_v(ii  ,jj  ,kk+1) id_cen(ii,jj,kk+1)];
        end
    end
end



end
function [cells, vertices] = convert_sgrid_to_prism_ugrid(xv,yv,zv_top,zv_vertical)

%
%  IDs of prismatic cells that are generated from structured grid
%
%          o ------- o ------- o
%          | \  19 / | \  23 / |
%          |  \   /  |  \   /  |
%          |20  o  18|24  o  22|
%          |  /   \  |  /   \  |
%          | /  17  \| /  21 \ |
%          o ------- o ------- o
%          | \  11 / | \ 15  / |
%          |  \   /  |  \   /  |
%          |12  o  10|17  o  14|
%          |  /   \  |  /   \  |
%          | /  9  \ | / 13  \ |
%          o ------- o ------- o
%          | \  3  / | \  7  / |
%          |  \   /  |  \   /  |
%          | 4  o  2 | 8  o  6 |
%   y      |  /   \  |  /   \  |
%  /|\     | /  1  \ | /  5  \ |
%   |      o ------- o ------- o
%   |
%    -------> x
%

[nvx,nvy] = size(xv);

nx = nvx-1;
ny = nvy-1;

nvz = length(zv_vertical);
nz  = nvz - 1;

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
            vertices(count,3) = zv_top(ii,jj) - zv_vertical(nvz-kk+1);
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
                zv_top(ii  ,jj  ) + ...
                zv_top(ii+1,jj  ) + ...
                zv_top(ii+1,jj+1) + ...
                zv_top(ii  ,jj+1) )/4 - zv_vertical(nvz-kk+1);
            
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
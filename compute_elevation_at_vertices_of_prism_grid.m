function zv_top = compute_elevation_at_vertices_of_prism_grid(sgrid, zc_top)

nx = sgrid.nx;
ny = sgrid.ny;

zv_top = zeros(nx+1,ny+1);

zv_top(2:nx,2:ny) = ...
    (...
    zc_top(1:nx-1,1:ny-1) + ...
    zc_top(2:nx  ,1:ny-1) + ...
    zc_top(1:nx-1,2:ny  ) + ...
    zc_top(2:nx  ,2:ny  ) ...
    )/4;

zv_top(2:nx,1   ) = (zc_top(1:nx-1,1     ) + zc_top(2:nx  ,1   ))/2;
zv_top(2:nx,ny+1) = (zc_top(1:nx-1,ny    ) + zc_top(2:nx  ,ny  ))/2;
zv_top(1   ,2:ny) = (zc_top(1     ,1:ny-1) + zc_top(1     ,2:ny))/2;
zv_top(nx+1,2:ny) = (zc_top(nx    ,1:ny-1) + zc_top(nx    ,2:ny))/2;

zv_top(1   ,1   ) = zc_top(1 ,1 );
zv_top(1   ,ny+1) = zc_top(1 ,ny);
zv_top(nx+1,1   ) = zc_top(nx,1 );
zv_top(nx+1,ny+1) = zc_top(nx,ny);

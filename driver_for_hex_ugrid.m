function driver_for_hex_ugrid()

sgrid.nx = 12;
sgrid.ny = 12;
sgrid.nz = 525;
sgrid.dx = 5000;
sgrid.dy = 5000;
sgrid.dz = 2;
sgrid.origin_x = 0;
sgrid.origin_y = 0;
sgrid.origin_z = 50;

sgrid.h5_material_filename = '/Users/bish218/projects_1drive/elm-pflotran-coupling/dataset/5km/pflotran/model_bcs/HFR_material_river.h5';
sgrid.h5_region_filename = '/Users/bish218/projects_1drive/elm-pflotran-coupling/dataset/5km/pflotran/model_bcs/HFR_material_river.h5';
h5_ugrid_filename = '/Users/bish218/projects_1drive/elm-pflotran-coupling/dataset/5km/pflotran/5km_hex_ugrid.h5';

convert_sgrid_to_ugrid(sgrid,h5_ugrid_filename);


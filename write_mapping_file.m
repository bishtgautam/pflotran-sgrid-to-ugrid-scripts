function write_mapping_file(filename, elm_nsoil, mapping)

disp(filename)
fid=fopen(filename,'w');
fprintf(fid,'! Num of CLM soil layers (for soil hydrology)\n');
fprintf(fid,'clm_nlevsoi %d\n\n',elm_nsoil);
fprintf(fid,'! Num of CLM ground layers (for soil heat transport)\n');
fprintf(fid,'clm_nlevgrnd 15\n\n');
fprintf(fid,'! Num of CLM layers mapped\n');
fprintf(fid,'clm_nlev_mapped %d\n\n',elm_nsoil);
fprintf(fid,'! Num of PFLOTRAN soil layers\n');
fprintf(fid,'pflotran_nlev %d\n\n',elm_nsoil);
fprintf(fid,'! Num of PFLOTRAN soil layers mapped\n');
fprintf(fid,'pflotran_nlev_mapped %d\n\n',elm_nsoil);
fprintf(fid,'! Num of weights\n');
fprintf(fid,'num_weights %d\n\n',size(mapping,1));
fprintf(fid,'! FORMAT:\n');
fprintf(fid,'! <DESTINATION_GRID>_cell_idx <SOURCE_GRID>_cell_idx weight\n');
fprintf(fid,'!\n');
fprintf(fid,'! Note: %s<DESTINATION_GRID>_cell_idx%s and %s<SOURCE_GRID>_cell_idx%s are in natural-order\n',char(39),char(39),char(39),char(39));

fprintf(fid,'!\n');for ii=1:size(mapping,1)
    fprintf(fid,'%d %d %f\n',mapping(ii,:));
end
fclose(fid);

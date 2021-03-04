function [K, dist_param, rows_cols] = load_parameters (path)

  file_lines = dlmread(path);
  K = file_lines(2:4,1:3);
  dist_param = file_lines(5,5:8)';
  rows_cols = file_lines(6:7,3);

endfunction

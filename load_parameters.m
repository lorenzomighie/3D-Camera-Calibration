function [K_truth, dist_param_truth] = load_parameters (path)

  ground_truth = dlmread(path);
  K_truth = ground_truth(2:4,1:3);
  dist_param_truth = ground_truth(5,5:8);
  rows_cols = ground_truth(6:7,3);

endfunction

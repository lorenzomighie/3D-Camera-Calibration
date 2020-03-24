function [K, d_param, n_inliers, chi_stats] = calibration_ls (K_g, d_param_g, w_p_3D, c_R_w, i_p_2D, iterations, damping, kernel_threshold)
  
  % example
  points = w_p_3D(1, :);
  points_h = [points'; 1];
  
  
  R = squeeze(c_R_w(1, :, :)); % R of exp 1
  z = squeeze(i_p_2D(1, 1, :)) % first point of exp 1
  proj_point = R*points_h;

  [e, J] = error_and_jacobian(K_g, d_param_g, proj_point, z);
 
  K = 0;
  d_param = 0;
  n_inliers = 0;
  chi_stats = 0; 
  
endfunction

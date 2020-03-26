function [K, d_param, n_inliers, chi_stats] = calibration_ls (K_g, d_param_g, w_p, c_R_w, i_p_2D, n_iterations, damping, kernel_threshold)
 
  K = zeros(3);
  d_param = zeros(1,4);
 
  n_p = length(w_p);
  state_dim = 8;
  
  chi_stats=zeros(1,n_iterations);
  n_inliers=zeros(1,n_iterations);
  p_hom = ones(n_p, 4);
  p_hom(:, 1:3) = w_p; % all world points set to homogeneous

  
  for it = 1:n_iterations
  
    H = zeros(state_dim, state_dim);
    b = zeros(state_dim, 1);
    
    for i_measurement = 1:length(i_p_2D)
      
      H = zeros(state_dim, state_dim);
      b = zeros(state_dim, 1);
      
      
      R = squeeze(c_R_w(i_measurement, :, :)); % R of exp m
            
      % function that iterates over all points
      [e, J, n_in, chi, H, b] = error_and_jacobian(p_hom, R, K_g, d_param_g, i_p_2D, i_measurement, n_p, state_dim, kernel_threshold, H, b);  
      n_inliers(it) = n_in;
      chi_stats(it) += chi;
      
      
      H+=eye(8)*damping;
      dx = -H\b; % 8 x 1
    
      [fx, fy, cx, cy, k1, k2, p1, p2] = get_X(K, d_param_g);
    
      % update
      K_g(1,1) = fx + dx(1);
      K_g(2,2) = fy + dx(2);
      K_g(1,3) = cx + dx(3);
      K_g(2,3) = cy + dx(4);
      d_param_g(1) = k1 + dx(5);
      d_param_g(2) = k2 + dx(6);
      d_param_g(3) = p1 + dx(7);
      d_param_g(4) = p2 + dx(8);
      
      
    endfor

    
    
  endfor
 
  K = K_g;
  d_param = d_param_g;
   
endfunction

function [K, d_param, n_inliers, chi_stats] = calibration_ls (K_g, d_param_g, w_p, c_R_w, i_p_2D, n_iterations, damping, kernel_threshold)
 
  n_inliers = 0;
  chi_stats = 0; 
  
  n_p = length(w_p);
  p_hom = ones(n_p, 4);
  p_hom(:, 1:3) = w_p; % all world points set to homogeneous
  
  R = squeeze(c_R_w(1, :, :)); % R of exp 1
  
  state_dim = 8;
  H = zeros(state_dim, state_dim);
  b=zeros(state_dim, 1);
  J = zeros(2*n_p, state_dim);
  e = zeros(2*state_dim, 1);
  chi_stats=zeros(1,n_iterations);
  n_inliers=zeros(1,n_iterations);
  
  for it = 1:n_iterations
  
    for i_measurement = 1:length(i_p_2D)
      
      for i_point = 1:length(w_p)
        
        z = squeeze(i_p_2D(1, i_point, :)); % point i of exp 1
        cam_p = (R*p_hom(i_point, :)')(1:3);
        
        [e, J_point] = error_and_jacobian(K_g, d_param_g, cam_p, z); 
        
        % chi error with threshold
        chi=e'*e;
        if (chi>kernel_threshold)
          e*=sqrt(kernel_threshold/chi);
          chi=kernel_threshold;
        else
          n_inliers(it)++;
        endif
        
        % fill the error 
        e(2*i_point - 1 :2*i_point) = e;
          
        % fill the jacobian
        J(2*i_point - 1 :2*i_point, :) = J_point;
      
      endfor
      
      chi_stats(it) += chi;
      H += J'*J;
      b += J'*e;
      
    endfor

    % H+=eye(6)*damping;
    dx = -H\b; % 8 x 1
    
    fx = K_g(1,1);
    fy = K_g(2,2);
    cx = K_g(1,3);
    cy = K_g(2,3);
    k1 = d_param_g(1);
    k2 = d_param_g(2);
    p1 = d_param_g(3);
    p2 = d_param_g(4);
    
    K_g(1,1) = fx + dx(1);
    K_g(2,2) = fy + dx(2);
    K_g(1,3) = cx + dx(3);
    K_g(2,3) = cy + dx(4);
    d_param_g(1) = k1 + dx(5);
    d_param_g(2) = k2 + dx(6);
    d_param_g(3) = p1 + dx(7);
    d_param_g(4) = p2 + dx(8);
    %X=v2t(dx)*X;
 
  endfor
 
  K = K_g;
  d_param = d_param_g;
  
  
endfunction

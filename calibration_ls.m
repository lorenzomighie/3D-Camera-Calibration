function [K, d_p, n_inliers, chi_stats] = calibration_ls (K_g, d_p_g, w_p, c_R_w, i_p, n_it, damping, k_t)
 
  % initialization
  K = K_g;
  d_p = d_p_g;
  n_p = length(w_p);
  state_dim = 8;
  chi_stats=zeros(1,n_it);
  n_inliers=zeros(1,n_it);
  % set all world points to homogeneous
  w_p_hom = ones(n_p, 4);
  w_p_hom(:, 1:3) = w_p; 

  for it = 1:n_it
  
    H = zeros(state_dim, state_dim);
    b = zeros(state_dim, 1);
    
    %for i_m = 1:length(i_p)
    for i_m = 1:length(i_p)
      
      % [R,t] of measurement m
      R_t = squeeze(c_R_w(i_m, :, :));  
      i_p_m = squeeze(i_p(i_m, :, :)); 
             
      % function that iterates over all points
      [e, J] = error_and_jacobian(w_p_hom, R_t, K, d_p, i_p_m, n_p, state_dim);  
       
      chi = e'*e;
      if (chi>k_t)
        e *= sqrt(k_t/chi);
        chi = k_t;
      else
        n_inliers(it)++;
      endif
      chi_stats(it) += chi;
      H += J'*J;
      b += J'*e;
      %H+=eye(8)*damping;
      
    endfor
    
    dx = -H\b; % 8 x 1
    
    [fx, fy, cx, cy, k1, k2, p1, p2] = get_X(K, d_p);
    
    % update
    K(1,1) = fx + dx(1);
    K(2,2) = fy + dx(2);
    K(1,3) = cx + dx(3);
    K(2,3) = cy + dx(4);
    d_p(1) = k1 + dx(5);
    d_p(2) = k2 + dx(6);
    d_p(3) = p1 + dx(7);
    d_p(4) = p2 + dx(8);
    
  endfor
   
endfunction

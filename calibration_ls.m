function [K, d_param, n_inliers, chi_stats] = calibration_ls (K_g, d_param_g, w_p, c_R_w, i_p_2D, iterations, damping, kernel_threshold)
  
  
  
  % for 1:iteration
  %   for 1:n_measurements
  %     [e,J] = errorAndJacobianManifold(X, P(:,i), Z(:,i));
  %     chi=e'*e;
  %     if (chi>kernel_threshold)
	%       e*=sqrt(kernel_threshold/chi);
	%       chi=kernel_threshold;
  %     else
	%       num_inliers(iteration)++;
  %     endif;
  %     chi_stats(iteration)+=chi;
  %     H+=J'*J;
  %     b+=J'*e;
  %   endfor
  %   H+=eye(6)*damping;
  %   dx=-H\b;
  
  % example
  p_hom = ones(length(w_p), 4);
  p_hom(:, 1:3) = w_p; % all world points set to homogeneous
  
  R = squeeze(c_R_w(1, :, :)) % R of exp 1
  z = squeeze(i_p_2D(1, :, :)); % all points of exp 1
  
  cam_p = zeros(length(w_p), 4);
  
  for index = 1:length(w_p)
    cam_p(index, :) = R*p_hom(index, :)'; % get all points in camera frame
  endfor
  cam_p
  %[e, J] = error_and_jacobian(K_g, d_param_g, proj_point, z);
 
  K = 0;
  d_param = 0;
  n_inliers = 0;
  chi_stats = 0; 
  
endfunction

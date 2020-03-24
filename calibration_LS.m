function [K, d_param, n_inliers, chi_stats] = calibration_LS (K_g, d_param_g, w_p_3D, c_R_w, i_p_2D, iterations, damping, kernel_threshold)
  
  function [e, J] = error_and_Jacobian(K, d_p, x_projected, z) % x_proj is c_R_w*[w_p_3D,1]
   % BEWARE: SHOULD WORK WITH VECTORS (entire set of points for a measurement)
   xi = x_projected(1)/x_projected(3);
   yi = x_projected(2)/x_projected(3);
   fx = K(1,1);
   fy = K(2,2);
   cx = K(1,3);
   cy = K(2,3);
   [k1, k2, p1, p2] = d_param_g;
   r2 = xi^2 + yi^2;
   xii = xi*(1 + k1*r2 + k2*r2^2) + 2*p1*xi*yi + p2*(r2 + 2*xi^2);
   yii = yi*(1 + k1*r2 + k2*r2^2) + p1*(r2 + 2*yi^2) + 2*p2*xi*yi;
   
   % prediction
   u = fx * xii + cx;
   v = fy * yuu + cy;
   
   e = zeros(2);
   e = [u, v] - z;
   
   J = zeros(2, 8);
   J(1, 1) = xii;
   J(2, 2) = yii;
   J(3:4, 3:4) = eye(2);
   J(1, 5:6) = fx*xi*r2*[1, r2];
   J(2, 5:6) = fy*yi*r2*[1, r2];
   J(1, 7) = fx*2*xi*yi;
   J(1, 8) = fx*(r2 + 2*xi^2);
   J(2, 7) = fy*(r2 + 2*yi^2);
   J(2, 8) = fy*xi*yi;
   
  endfunction

endfunction

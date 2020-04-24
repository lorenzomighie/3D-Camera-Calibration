function [e, J] = error_and_jacobian(p_hom, R, K, d_p, i_p_m, n_p, s_dim, r_c) 

  J = zeros(2*n_p, s_dim);
  e = zeros(2*n_p, 1);
  
  % for taking into account point out of bounds
  index_skip = []; 
  rows = r_c(1);  % 640 
  cols = r_c(2);  % 480
  
  for i_point = 1:n_p
    z = squeeze(i_p_m(i_point, :));
    cam_p = (R*p_hom(i_point, :)')(1:3);  
    
    % projection
    xi = cam_p(1)/cam_p(3);
    yi = cam_p(2)/cam_p(3);
    [fx, fy, cx, cy, k1, k2, p1, p2] = get_X(K, d_p);
    
    % measurement function with Brown–Conrady distortion model
    r2 = xi^2 + yi^2;
    xii = xi*(1 + k1*r2 + k2*r2^2) + 2*p1*xi*yi + p2*(r2 + 2*xi^2);
    yii = yi*(1 + k1*r2 + k2*r2^2) + p1*(r2 + 2*yi^2) + 2*p2*xi*yi;
    u = fx * xii + cx;
    v = fy * yii + cy;
   
    % check image limits
    if(u < 0 || u > cols || v < 0 || v > rows)
      skip = 'point out of bound'
      % track index of points out of bounds
      index_skip = [index_skip, i_point]
    endif
   
    % point error
    e_p = [u; v] - z';
    % point jacobian
    J_p = zeros(2, 8);
    J_p(1, 1) = xii;
    J_p(2, 2) = yii;
    J_p(1:2, 3:4) = eye(2);
    J_p(1, 5:6) = fx*xi*r2*[1, r2];
    J_p(2, 5:6) = fy*yi*r2*[1, r2];
    J_p(1, 7) = fx*2*xi*yi;
    J_p(1, 8) = fx*(r2 + 2*xi^2);
    J_p(2, 7) = fy*(r2 + 2*yi^2);
    J_p(2, 8) = fy*2*xi*yi;
        
    % fill measurement error and jacobian
    e(2*i_point - 1 :2*i_point) = e_p;   
    J(2*i_point - 1 :2*i_point, :) = J_p;
    
  endfor
  
  % remove point jacobian and error out of boundS
  for i_skip = 1:length(index_skip)
    i = index_skip(i_skip)
    e(2*i-1:2*i) = [];
    J(2*i - 1 :2*i, :) = [];
    index_skip -= 1;
    
  endfor 
endfunction





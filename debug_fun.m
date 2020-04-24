

function [retval] = debug_fun (measurement, point, w_p_3D, c_R_w, K_t, d_param_t,
                               i_p_2D)

  point_n = w_p_3D(point,:)' ;
  R = squeeze(c_R_w(measurement,:,:)) ;
  point_n_c = R*[point_n; 1]; 
  xi = point_n_c(1)/point_n_c(3);
  yi = point_n_c(2)/point_n_c(3); 
  r2 = xi^2 + yi^2;
  [fx, fy, cx, cy, k1, k2, p1, p2] = get_X(K_t, d_param_t); 
  xii = xi + xi*( k1*r2 + k2*r2^2) + 2*p1*xi*yi + p2*(r2 + 2*(xi)^2);
  yii = yi + yi*( k1*r2 + k2*r2^2) + p1*(r2 + 2*(yi)^2) + 2*p2*xi*yi;
  u = fx * xii + cx;
  v = fy * yii + cy;
  point_n_ima_truth = [u; v]
  observation_n = squeeze(i_p_2D(measurement,point,:))

endfunction

function [e, J] = error_and_jacobian(K, d_p, cam_p, z) 
  
   xi = cam_p(1)/cam_p(3);
   yi = cam_p(2)/cam_p(3);
   fx = K(1,1);
   fy = K(2,2);
   cx = K(1,3);
   cy = K(2,3);
   k1 = d_p(1);
   k2 = d_p(2);
   p1 = d_p(3);
   p2 = d_p(4);
   r2 = xi^2 + yi^2;
   xii = xi*(1 + k1*r2 + k2*r2^2) + 2*p1*xi*yi + p2*(r2 + 2*xi^2);
   yii = yi*(1 + k1*r2 + k2*r2^2) + p1*(r2 + 2*yi^2) + 2*p2*xi*yi;
   
   % prediction
   u = fx * xii + cx;
   v = fy * yii + cy;
   
   e = [u; v] - z;
   
   J = zeros(2, 8);
   J(1, 1) = xii;
   J(2, 2) = yii;
   J(1:2, 3:4) = eye(2);
   J(1, 5:6) = fx*xi*r2*[1, r2];
   J(2, 5:6) = fy*yi*r2*[1, r2];
   J(1, 7) = fx*2*xi*yi;
   J(1, 8) = fx*(r2 + 2*xi^2);
   J(2, 7) = fy*(r2 + 2*yi^2);
   J(2, 8) = fy*xi*yi;
   
endfunction

function [e, J, skip_point] = error_and_jacobian(K, d_p, cam_p, z) 

   xi = cam_p(1)/cam_p(3);
   yi = cam_p(2)/cam_p(3);
   ima_p = K*cam_p;
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
   
   % TEMPORARY: TRY WITHOUT DISTORTION COEFF
   %u = fx * xi + cx;
   %v = fy * yi + cy;
   % TEMPORARY: TRY WITHOUT DISTORTION COEFF
   
   
   if(u < 0 || u > 640 || v < 0 || v > 480)
    skip_point = true;
   else
    skip_point = false;
   endif
   
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
   J(2, 8) = fy*2*xi*yi;
   
   % TEMPORARY: TRY WITHOUT DISTORTION COEFF
   %a = [fx*cam_p(1) + cx*cam_p(3); fy*cam_p(2) + cy*cam_p(3); cam_p(3)];
   %dproj_da = [1/a(3), 0, -a(1)/a(3)^2;
   %            0, 1/a(3), -a(2)/a(3)^2];
   %da_dx = [cam_p(1), 0, cam_p(3), 0;
   %         0, cam_p(2), 0, cam_p(3);
   %         0,    0,     0,   0];
   %         
   %Jt = dproj_da*da_dx;
   
   %Jt = [ xi, 0, 1, 0;
   %       0, yi, 0, 1];
  
   % TEMPORARY: TRY WITHOUT DISTORTION COEFF
   
endfunction

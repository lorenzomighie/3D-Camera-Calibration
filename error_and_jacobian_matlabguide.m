function [e, J, skip_point] = error_and_jacobian_matlabguide(K, d_p, cam_p, z) 

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
   
   
   % prediction
   u_nd = fx * xi + cx;
   v_nd = fy * yi + cy;
   
   u_d = u_nd*(1 + k1*r2 + k2*r2^2) + 2*p1*u_nd*v_nd + p2*(r2 + 2*u_nd^2);
   v_d = v_nd*(1 + k1*r2 + k2*r2^2) + p1*(r2 + 2*v_nd^2) + 2*p2*u_nd*v_nd;
   
   % TEMPORARY: TRY WITHOUT DISTORTION COEFF
   %u = fx * xi + cx;
   %v = fy * yi + cy;
   % TEMPORARY: TRY WITHOUT DISTORTION COEFF
   
   
   if(u_d < 0 || u_d > 640 || v_d < 0 || v_d > 480)
    skip_point = true;
   else
    skip_point = false;
   endif
   
   e = [u_d; v_d] - z;
   
   J = zeros(2, 8);
   J(1, 1) = u_nd*(1+k1*r2+k2*r2^2) + 2*p1*u_nd*(fy*v_nd+cy) + 4*p2*u_nd*(fx*u_nd+cx);
   J(2, 2) = v_nd*(1+k1*r2+k2*r2^2) + 2*p2*v_nd*(fx*u_nd+cx) + 4*p1*v_nd*(fy*v_nd+cy);
   J(1 ,2) = 2*p1*v_nd*(fx*u_nd+cx);
   J(2 ,1) = 2*p2*u_nd*(fy*v_nd+cy);
   J(1, 3) = (1+k1*r2+k2*r2^4) + 2*p1*(fy*v_nd+cy) + 4*p2*(fx*u_nd+cx);
   J(2, 4) = (1+k1*r2+k2*r2^4) + 2*p2*(fx*u_nd+cx) + 4*p1*(v_nd*fy+cy);
   J(1, 4) = 2*p1*(fx*u_nd+cx); 
   J(2, 3) = 2*p2*(fy*v_nd+cy);
   J(1, 5) = r2*(fx*u_nd+cx);
   J(2, 5) = r2*(fy*v_nd+cy);
   J(1, 6) = r2^2*(fx*u_nd+cx);
   J(2, 6) = r2^2*(fy*v_nd+cy);
   J(1, 7) = 2*(fx*u_nd+cx)*(fy*v_nd+cy);
   J(2, 8) = 2*(fx*u_nd+cx)*(fy*v_nd+cy);
   J(1, 8) = r2 + 2*(u_nd*fx+cx)^2;
   J(2, 7) = r2 + 2*(v_nd*fy+cy)^2;
   

   
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
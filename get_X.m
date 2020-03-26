function [fx, fy, cx, cy, k1, k2, p1, p2] = get_X(K, d_p)

  fx = K(1,1);
  fy = K(2,2);
  cx = K(1,3);
  cy = K(2,3);
  k1 = d_p(1);
  k2 = d_p(2);
  p1 = d_p(3);
  p2 = d_p(4);
  
endfunction
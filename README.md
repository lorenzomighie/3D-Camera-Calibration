# 3D-Camera-Calibration

Calibration of a pinhole camera with Brown-Conrady distortion model using Gauss-Newton optimization 

# Background

## Pinhole Camera Model:

K = [fx, 0, cx;
     0, fy, cy;
     0,  0,  1]
     
[x_im; y_im] = K*[x_cam; y_cam; z_cam]
 
In which the vector X_cam = [x_cam; y_cam; z_cam] is expressed with respect to the camera coordinate system, and X_im = [x_im; y_im] contains the points 
coordinates in the image plane.

Furthermore

X_cam' = w_R_c*X_w' 

In which the vector ' are in homogeneous coordinates and the matrix w_R_c is the rota-translation matrix expressing the pose of the camera with regard to 
the world frame.
 
## Brown-Conrady distortion model

It's a simple distortion model that takes into account radial and tangential distortion with 4 parameters (k1, k2, p1 and p2); the following formula show how 
to find the image coordinates with it.

xi = x_cam(1)/y_cam(3);
yi = x_cam(2)/y_cam(3);

r2 = xi^2 + yi^2;
xii = xi*(1 + k1*r2 + k2*r2^2) + 2*p1*xi*yi + p2*(r2 + 2*xi^2);
yii = yi*(1 + k1*r2 + k2*r2^2) + p1*(r2 + 2*yi^2) + 2*p2*xi*yi;

u = fx * xii + cx;
v = fy * yii + cy;

## State Definition

Such an optimization task requires non-linear least squares optimization applied to 8 parameters: fx, fy, cx, cy, p1, p2, d1, d2, therefore the state vector
is the vector containing this 8 elements.

## Data available

For this calibration we have available: 
   - the 3D point position in world frame (arranged on a chessboard)
   - the initial state parameters (their initial guess)
   - a set of measurements, each one containing:
      * the pose of the camera with respect to the world;
      * the points measured.
      
The association of the 3D point with the respective measurements is known.

## Gauss Newton Procedure


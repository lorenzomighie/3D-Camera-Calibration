# 3D-Camera-Calibration

Calibration of a pinhole camera with Brown-Conrady distortion model using Gauss-Newton optimization with Octave

# Background

## Pinhole Camera Model:

K = [fx, 0, cx;
     0, fy, cy;
     0,  0,  1].
     
[x_im; y_im] = K [x_cam; y_cam; z_cam].
 
In which the vector X_cam = [x_cam; y_cam; z_cam] is expressed with respect to the camera coordinate system, and X_im = [x_im; y_im] contains the points 
coordinates in the image plane.

Furthermore:

X_cam' = w_R_c X_w' 

In which the vector ' are in homogeneous coordinates and the matrix w_R_c is the rota-translation matrix expressing the pose of the camera with regard to 
the world frame.
 
## Brown-Conrady distortion model

It's a simple distortion model that takes into account radial and tangential distortion with 4 parameters (k1, k2, p1 and p2); the following formula show how 
to find the image coordinates with it.

xi = x_cam/z_cam;

yi = y_cam/z_cam;

r2 = xi^2 + yi^2;

xii = xi (1 + k1 r2 + k2 r2^2) + 2 p1 xi yi + p2*(r2 + 2 xi^2);

yii = yi (1 + k1 r2 + k2 r2^2) + p1 (r2 + 2*yi^2) + 2 p2 xi yi;

u = fx xii + cx;

v = fy yii + cy.

# State Definition

Such an optimization task requires non-linear least squares optimization applied to 8 parameters: fx, fy, cx, cy, p1, p2, d1, d2, therefore the state vector X
is the vector containing this 8 elements.

# Data available

For this calibration we have available: 
   - the 3D point position in world frame (arranged on a chessboard);
   - the initial state parameters (their initial guess);
   - a set of measurements, each one containing:
      * the pose of the camera with respect to the world;
      * the points measured.
      
In addition, the association of the 3D point with the respective measurements is known.

# Gauss Newton Procedure

To apply this algorithm it is required to:

   - compute the error for each measurement m, given as the sum of the error of each point measured: e_i = (u_i, v_i) - (z_x_i, z_y_i), where the vector 
Z_i is the i-th point of the measurement m;
   - compute the Jacobian for each measurement m, by deriving the error with respect to the state, for each point: J_i = d(e_i)/dX.

For each measurement an error vector e of size (2 number_of_points) and a Jacobian J of size (2 number_of_points, space_dimension = 8) are computed, and the matrices H and the vector b are updated: H += J'*J; b += J'*e.

Once this procedure has been done for each measurement, the optimization step takes places: dx = -H\b --> X = X + dx. All the aforementioned steps can be done for  a number of iterations, after which the algorithm will converge.

# How To Use

Just run the file 'main.m' with Octave.

If you want to use different camera model you'll have to modify:
- the state vector;
- the measurement function, hence the error vector;
- the Jacobian;
- the size of the matrix H and the vector b.

# Results

The algorithm reaches convergence after 4 iterations, with an Eucledean norm of about 0.031.

![GitHub Logo](/images/error_plot.png)
Format: ![Alt Text](url)


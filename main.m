%% clear history

clc
clear all
close all

%% set format of data to be dislayed as long

format long g

%% load dataset

data_dir = 'data/';
truth_file = 'camera_param_ground_truth.dat';
guess_file = 'camera_param_initial_guess.dat';
points_file = 'grid3d.dat';
measurements_initial_name = 'image';

[K_t, d_param_t, row_col] = load_parameters(strcat(data_dir, truth_file));
[K_g, d_param_g, ~] = load_parameters(strcat(data_dir, guess_file));

w_p_3D = dlmread(strcat(data_dir, points_file)); % points in world frame
n_points = size(w_p_3D)(1);

[c_R_w, i_p_2D, n_measurements] = load_measurements(data_dir, measurements_initial_name, n_points);
% set of world to camera Pose matrix and image points

% (at the end I could try with linear relaxation)

%% Run Least Squares
iterations=10;
damping=0; # damping factor
kernel_threshold = 1e9;

%[K_LS, d_param_LS,  n_inliers, error_chi] = calibration_ls(K_g, d_param_g, w_p_3D, c_R_w, i_p_2D, iterations, damping, kernel_threshold);

% notice that the 3d points give are on a plane and belong to a grid

% DEBUG --> see if observation corresponds to w_point projected with ground truth
measurement = 1
point = 1
point_n = w_p_3D(point,:)' 
R = squeeze(c_R_w(measurement,:,:)) 
K_g

point_n_c = R*[point_n; 1]; 
point_n_im = K_g*point_n_c(1:3); 
point_n_ima_guess = point_n_im(1:2)/point_n_im(3)


[fx, fy, cx, cy, k1, k2, p1, p2] = get_X(K_t, d_param_t); %g
xi = point_n_c(1)/point_n_c(3);
yi = point_n_c(2)/point_n_c(3);
K_t
d_param_t
r2 = xi^2 + yi^2;
xii = xi*(1 + k1*r2 + k2*r2^2) + 2*p1*xi*yi + p2*(r2 + 2*xi^2);
yii = yi*(1 + k1*r2 + k2*r2^2) + p1*(r2 + 2*yi^2) + 2*p2*xi*yi;
u = fx * xii + cx;
v = fy * yii + cy;
point_n_ima_truth = [u; v]

observation_n = squeeze(i_p_2D(measurement,point,:))

%% plot error
%figure
%plot([1:iterations], error_chi)
%K_LS
%K_t
%d_param_LS
%d_param_t
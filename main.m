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
% set of Pose matrix world to camera and image points

%% should add a check whether points are inside the rows and columns wrt to the offset
% also check if camera is pointing in the right direction
% (at the end I could try with linear relaxation)
% Using the model with radial and tangential distortion
%% Run Least Squares

iterations=100;
damping=0; # damping factor
kernel_threshold = 1e9;

[K_LS, d_param_LS,  n_inliers, error_chi] = calibration_ls(K_t, d_param_t, w_p_3D, c_R_w, i_p_2D, iterations, damping, kernel_threshold);

% notice that the 3d points give are on a plane and belong to a grid

%% plot error
figure
plot([1:iterations], error_chi)
K_LS
K_t
d_param_LS
d_param_t
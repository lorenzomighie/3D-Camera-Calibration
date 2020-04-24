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

% set of world to camera Pose matrix and image points
[c_R_w, i_p_2D, n_measurements] = load_measurements(data_dir, measurements_initial_name, n_points);

%% Run Least Squares
iterations=20;
kernel_threshold = 1e9; 
state_dim = 8;

[K_LS, d_param_LS,  n_inliers, error_chi] = calibration_ls(state_dim, K_g, d_param_g, 
  w_p_3D, c_R_w, i_p_2D, iterations, kernel_threshold, row_col);

% print parameters obtained
K_LS
K_t
d_param_LS
d_param_t

% plot error
figure
plot([1:iterations], error_chi)
final_error = error_chi(end)

% DEBUG --> see if observation corresponds to w_point projected with ground truth
%measurement = 3
%point = 7
%debug_fun(measurement, point, w_p_3D, c_R_w, K_t, d_param_t, i_p_2D)
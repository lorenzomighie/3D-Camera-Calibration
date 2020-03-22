%% clear history

clc
clear all
close all

%% load dataset

data_dir = 'data/';
truth_file = 'camera_param_ground_truth.dat';
guess_file = 'camera_param_initial_guess.dat';
points_file = 'grid3d.dat';
measurements_initial_name = 'image';

[K_t, d_param_t] = load_parameters(strcat(data_dir, truth_file));
[K_g, d_param_g] = load_parameters(strcat(data_dir, guess_file));

w_p_3D = dlmread(strcat(data_dir, points_file)); % points in world frame
n_points = size(w_p_3D)(1)

[cRw, i_p_2D, n_measurements] = load_measurements(data_dir, measurements_initial_name, n_points);
% set of Pose matrix world to camera and image points

%% PROBLEM: THE NUMBER STORED HAVE ONLY 4 CIFRE !!!! --> use double, recheck function 
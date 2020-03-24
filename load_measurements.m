function [camera_R_world, image_p_2D, n_measurements] = load_measurements (data_dir, name, n_points)

  n_measurements = 0;
  file_list = readdir(data_dir);
  n_el = numel(file_list);
  camera_R_world = zeros(n_el, 4, 4); % at most, the first dimension will be truncated later
  image_p_2D = zeros(n_el, n_points, 2);

  for n_file=3:n_el
    file = file_list{n_file};
    if (file(1:5) == name)
      n_measurements += 1;
      path_file = strcat(data_dir, file);
      read_file = dlmread(path_file);
      camera_R_world(n_measurements, :, :) = read_file(2:5, 1:4);
      image_p_2D(n_measurements, :, :) = read_file(7:end, 1:2); 
    endif
  endfor
  
  camera_R_world = camera_R_world(1:n_measurements, :, :);
  image_p_2D = image_p_2D(1:n_measurements, :, :);
  

endfunction

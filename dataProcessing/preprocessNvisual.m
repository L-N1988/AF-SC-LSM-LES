clc; clear; close all;
% Define base directory
base_dir = './velocityPlaneX0';

% Output folder
output_folder = 'centerPlaneU';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Get list of subdirectories
dir_list = dir(base_dir);
time_dirs = {dir_list([dir_list.isdir]).name}';
time_dirs = time_dirs(~ismember(time_dirs, {'.', '..'}));

% Sort numerically by converting to double
numeric_time_dirs = cellfun(@(x) str2double(x), time_dirs);
[~, idx] = sort(numeric_time_dirs);
time_dirs = time_dirs(idx);
% Store all data in a 3D array
time_dir = time_dirs{1};
fprintf('processing time = %s\n', time_dir);

% construct path to raw file
raw_path = fullfile(base_dir, time_dir, 'u_x0plane.raw');
if ~exist(raw_path, 'file')
    error('skipping %s: file not found - %s\n', time_dir, raw_path);
end

% read the raw file
fid = fopen(raw_path, 'r');
% skip the first line (metadata)
fgetl(fid);
% read second line (headers)
headers = fgetl(fid);
headers = strsplit(headers, ' ');

% read numerical data
data_lines = textscan(fid, '%f %f %f %f %f %f');
fclose(fid);

% convert to matrix
data = cell2mat(data_lines);
all_data = zeros([size(data), length(time_dirs)]);
all_data(:, :, 1) = data;

% Loop through each time directory
for i = 2:length(time_dirs)
    time_dir = time_dirs{i};
    fprintf('processing time = %s\n', time_dir);
    
    % construct path to raw file
    raw_path = fullfile(base_dir, time_dir, 'u_x0plane.raw');
    if ~exist(raw_path, 'file')
        fprintf('skipping %s: file not found - %s\n', time_dir, raw_path);
        continue;
    end
    
    % read the raw file
    fid = fopen(raw_path, 'r');
    % skip the first line (metadata)
    fgetl(fid);
    % read second line (headers)
    headers = fgetl(fid);
    headers = strsplit(headers, ' ');
    
    % read numerical data
    data_lines = textscan(fid, '%f %f %f %f %f %f');
    fclose(fid);
    
    % convert to matrix
    data = cell2mat(data_lines);
    all_data(:, :, i) = data;
    
    % Plotting function (nested inside loop)
    if mod(str2double(time_dir), 5) == 0
        % % Extract columns
        % x = data(:, 1);
        % y = data(:, 2);
        % z = data(:, 3);
        % Ux = data(:, 4);
        % Uy = data(:, 5);
        % Uz = data(:, 6);
        % 
        % % Create interpolation mesh
        % y_unique = unique(y);
        % z_unique = unique(z);
        % 
        % [Yq, Zq] = meshgrid(linspace(min(y), max(y), length(y_unique)*10), ...
        %     linspace(min(z), max(z), length(z_unique))*20);
        % 
        % % Interpolate using scatteredInterpolant (equivalent to CloughTocher2DInterpolator)
        % F_Ux = scatteredInterpolant(y, z, Ux, 'natural');
        % F_Uy = scatteredInterpolant(y, z, Uy, 'natural');
        % F_Uz = scatteredInterpolant(y, z, Uz, 'natural');
        % 
        % Ux_grid = F_Ux(Yq, Zq);
        % Uy_grid = F_Uy(Yq, Zq);
        % Uz_grid = F_Uz(Yq, Zq);
        % plot_and_save(Ux_grid, 'U_x', sprintf('Ux_%s.png', time_dir), Yq, Zq, output_folder);
        % plot_and_save(Uy_grid, 'U_y', sprintf('Uy_%s.png', time_dir), Yq, Zq, output_folder);
        % plot_and_save(Uz_grid, 'U_z', sprintf('Uz_%s.png', time_dir), Yq, Zq, output_folder);
    else
        continue;
    end
end

disp(['Plots saved in folder: ', output_folder]);
% Save all_data as a .mat file
save("all_data.mat", 'all_data', 'numeric_time_dirs');


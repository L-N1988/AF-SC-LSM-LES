clc; clear; close all;
data = load("all_data.mat");

Ts = data.numeric_time_dirs;
Fields = data.all_data;
% all sample data are in same plane
Xs = squeeze(Fields(:, 1, :));
Ys = squeeze(Fields(:, 2, :));
Zs = squeeze(Fields(:, 3, :));
x_unique = sort(unique(Xs(:, :)));
y_unique = sort(unique(Ys(:, :)));
z_unique = sort(unique(Zs(:, :)));
center_x = x_unique(round(length(x_unique)/2), 1);
center_y = y_unique(round(length(y_unique)/2), 1);
center_z = z_unique(round(length(z_unique)/2), 1);
% [velocity at point, time]
us = squeeze(Fields(:, 4, :));
vs = squeeze(Fields(:, 5, :));
ws = squeeze(Fields(:, 6, :));

%% Interpolate velocity field to get structured grid data

% Create interpolation mesh
[Zq, Yq] = meshgrid(linspace(min(z_unique), max(z_unique), length(z_unique)*10), ...
    linspace(min(y_unique), max(y_unique), length(y_unique)*20));

uq = zeros([size(Zq), length(Ts)]);
vq = zeros([size(Zq), length(Ts)]);
wq = zeros([size(Zq), length(Ts)]);

for i = 1:length(Ts)
    Fuq = scatteredInterpolant(Zs(:, 1), Ys(:, 1), us(:, i), 'linear');
    Fvq = scatteredInterpolant(Zs(:, 1), Ys(:, 1), vs(:, i), 'linear');
    Fwq = scatteredInterpolant(Zs(:, 1), Ys(:, 1), ws(:, i), 'linear');
    uq(:, :, i) = Fuq(Zq, Yq);
    vq(:, :, i) = Fvq(Zq, Yq);
    wq(:, :, i) = Fwq(Zq, Yq);
end

%% Velocity statistics
% Compute time-averaged velocity
U_t = mean(uq, 3);
V_t = mean(vq, 3);
W_t = mean(wq, 3);

% Compute double-averaged velocity (spatial average in spanwise z direction)
U_zt = mean(U_t, 2);
V_zt = mean(V_t, 2);
W_zt = mean(W_t, 2);

% Compute turbulent velocity components (fluctuations)
u_pri = uq - U_t; % u' = u - <u>
v_pri = vq - V_t; % v' = v - <v>
w_pri = wq - W_t; % w' = w - <w>

% Compute second moments
% stress tensor
uv = u_pri .* v_pri; % u'v'
uw = u_pri .* w_pri; % u'w'
vw = v_pri .* w_pri; % v'w'
% TKE components
uu = u_pri.^2; % u'u'
vv = v_pri.^2; % v'v'
ww = w_pri.^2; % w'w'

% Compute time-averaged second moments
uv_t  = mean(uv, 3);
uw_t  = mean(uw, 3);
vw_t  = mean(vw, 3);
RSS   = -1 * uv_t; % Reynolds shear stress

uu_t  = mean(uu, 3);
vv_t  = mean(vv, 3);
ww_t  = mean(ww, 3);
TKE = 0.5 * (uu_t + vv_t + ww_t); % Turbulent kinetic energy

% Compute root mean square (RMS) values for turbulence strength
u_rms = sqrt(uu_t); % RMS of u
v_rms = sqrt(vv_t); % RMS of v

% Compute double-averaged second moments (spatial average in spanwise z-direction)
uv_xt  = mean(uv_t, 2); % u'v'
uu_xt  = mean(uu_t, 2); % u'u'
vv_xt  = mean(vv_t, 2); % v'v'

%% Plot central line profiles
figure();
scatter(U_t(:, round(length(U_t)/ 2)), Yq(:, 1), 10, 'b', 'filled');
hold on;
scatter(V_t(:, round(length(V_t)/ 2)), Yq(:, 1), 10, 'r', 'filled');
scatter(W_t(:, round(length(W_t)/ 2)), Yq(:, 1), 10, 'g', 'filled');
xlabel('Velocity (m/s)');
ylabel('Y (m)');
legend('U', 'V', 'W');
title('Central line velocity profiles');
hold off;
saveas(gcf, 'uvw_central_line_profiles.png');

figure();
scatter(u_rms(:, round(length(u_rms)/ 2)), Yq(:, 1), 10, 'b', 'filled');
hold on;
scatter(v_rms(:, round(length(v_rms)/ 2)), Yq(:, 1), 10, 'r', 'filled');
scatter(w_rms(:, round(length(w_rms)/ 2)), Yq(:, 1), 10, 'g', 'filled');
xlabel('Turbulence intensity (m/s)');
ylabel('Y (m)');
legend('u_rms', 'v_rms');
title('Central line turbulence intensity profiles');
hold off;
saveas(gcf, 'rms_central_line_profiles.png');

figure();
scatter(uv_t(:, round(length(uv_t)/ 2)), Yq(:, 1), 10, 'b', 'filled');
hold on;
scatter(uw_t(:, round(length(uu_t)/ 2)), Yq(:, 1), 10, 'r', 'filled');
scatter(vw_t(:, round(length(vv_t)/ 2)), Yq(:, 1), 10, 'g', 'filled');
xlabel('Tensor (m^2/s^2)');
ylabel('Y (m)');
lengend('uv_xt', 'uw_xt', 'vw_xt');
title('Central line stress tensor profiles');
hold off;
saveas(gcf, 'tensor_central_line_profiles.png');

%% Plot 2D vector fields in Y-Z plane
plot_and_save(Ubar, '$\overline{U}$', 'Ux_2d.png', Yq, Zq, '.');
plot_and_save(Vbar, '$\overline{V}$', 'Uy_2d.png', Yq, Zq, '.');
plot_and_save(Wbar, '$\overline{W}$', 'Uz_2d.png', Yq, Zq, '.');
plot_and_save(u_rms, '$u_{rms}$', 'u_rms_2d.png', Yq, Zq, '.');
plot_and_save(v_rms, '$v_{rms}$', 'v_rms_2d.png', Yq, Zq, '.');
plot_and_save(w_rms, '$w_{rms}$', 'w_rms_2d.png', Yq, Zq, '.');
plot_and_save(RSS, '$RSS$', 'RSS_2d.png', Yq, Zq, '.');
plot_and_save(TKE, '$TKE$', 'TKE_2d.png', Yq, Zq, '.');

%% Plot V and W 2D vector fields in Y-Z plane
figur = figure;
contourf(Zq(1:10:end, 1:10:end), Yq(1:10:end, 1:10:end), sqrt(Vq(1:10:end, 1:10:end).^2 + Wq(1:10:end, 1:10:end).^2), 5);
hold on;
quiver(Zq(1:10:end, 1:10:end), Yq(1:10:end, 1:10:end), Wq(1:10:end, 1:10:end), Vq(1:10:end, 1:10:end), 'k', 'LineWidth', 1);
axis xy
colorbar;
colormap('sky')
set(gca, 'YDir', 'normal'); % Ensure correct orientation
xlabel('Y');
ylabel('Z');
title('2D Vector Field');
hold off;
saveas(figur, 'vector_field.png');

%% Save internal variables
% Save statistical results to files
Y = Yq(:, 1);
X = Zq(1, :);
xmesh = Zq;
ymesh = Yq;
save('u_stat.mat', 'xmesh', 'X', 'ymesh', 'Y', 'U_t', 'V_t', 'W_t', 'uv_t', 'uw_t', 'vw_t', 'u_rms', 'v_rms', 'w_rms', 'TKE');
save('u4pxx.mat', 'xmesh', 'X', 'ymesh', 'Y', 'U_t', 'V_t', 'u_pri', 'v_pri');
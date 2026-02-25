prop_data = readmatrix("12x6E.xlsx");

RPM = prop_data(:,1);
prop_vel = prop_data(:,2) /  2.237; % in m/s
prop_thrust = prop_data(:,12);

x = prop_vel(:);
y = RPM(:);
z = prop_thrust(:);

valid = isfinite(x) & isfinite(y) & isfinite(z);

x = x(valid);
y = y(valid);
z = z(valid);

Thrust_func = scatteredInterpolant(x, y, z);

[xq, yq] = meshgrid( ...
    linspace(min(x), max(x), 50), ...
    linspace(min(y), max(y), 50));

zq = Thrust_func(xq, yq);

figure()
scatter3(prop_vel, RPM, prop_thrust, 15, prop_thrust, 'filled')
hold on
xlabel("Velocity (m/s)")
ylabel("RPM")
set(gca, 'YDir', 'reverse')
zlabel("Thrust (N)")
colormap('parula')
colorbar

surf(xq, yq, zq, 'FaceAlpha', 0.5)
shading interp
colorbar
title("Multivariable Fit for Thrust vs. RPM and Velocity for APC 12x6E")
prop_dat = readmatrix("12x6E.xlsx");

RPM_np = prop_dat(:,1);
prop_vel_np = prop_dat(:,2) /  2.237; % in m/s
prop_np = prop_dat(:,4);

prop_vel_np = prop_vel_np(:);
prop_RPM_np = RPM_np(:);
prop_np = prop_np(:);

valid = isfinite(prop_vel_np) & isfinite(prop_RPM_np) & isfinite(prop_np);

prop_vel_np = prop_vel_np(valid);
prop_RPM_np = prop_RPM_np(valid);
prop_np = prop_np(valid);

etaP_func = scatteredInterpolant(prop_vel_np, prop_RPM_np, prop_np);

[vel_p, RPM_p] = meshgrid( ...
    linspace(min(prop_vel_np), max(prop_vel_np), 50), ...
    linspace(min(prop_RPM_np), max(prop_RPM_np), 50));

zq = etaP_func(vel_p, RPM_p);

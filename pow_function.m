prop_dats = readmatrix("12x6E.xlsx");

RPM_pow = prop_dats(:,1);
prop_vel_pow = prop_dats(:,2) /  2.237; % in m/s
prop_pow = prop_dats(:,10); % power in watts

prop_vel_pow = prop_vel_pow(:);
prop_RPM_pow = RPM_pow(:);
prop_pow = prop_pow(:);

valid = isfinite(prop_vel_pow) & isfinite(prop_RPM_pow) & isfinite(prop_pow);

prop_vel_pow = prop_vel_pow(valid);
prop_RPM_pow = prop_RPM_pow(valid);
prop_pow = prop_pow(valid);

pow_func = scatteredInterpolant(prop_vel_pow, prop_RPM_pow, prop_pow);

[vel_pow, RPM_pow] = meshgrid( ...
    linspace(min(prop_vel_pow), max(prop_vel_pow), 50), ...
    linspace(min(prop_RPM_pow), max(prop_RPM_pow), 50));

zq = etaP_func(vel_pow, RPM_pow);
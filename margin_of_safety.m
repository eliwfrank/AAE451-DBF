% margin of safety vs number of ribs plot

I = 1.2e-8;        % m^4  (example — replace)
y_max = 0.04;      % m    (distance to outer skin)

sigma_bend = M_bend .* y_max ./ I;   % 1x1000 vector

Q = 3.5e-6;       % m^3 (example)
t_web = 0.003;    % m

tau_web = V .* Q ./ (I * t_web);

nu = 0.3;
E_skin = 2*G_skin*(1+nu);
E_web  = 2*G_web*(1+nu);

N_ribs = 2:12;
semi_span = 0.9;    % example

k = 4;

MS_stringer_buckling = zeros(size(N_ribs));
MS_web_shear = zeros(size(N_ribs));
MS_bending_material = zeros(size(N_ribs));

sigma_allow = 120e6;     % define from material
tau_allow   = 20e6;

for i = 1:length(N_ribs)

    N = N_ribs(i);
    panel_length = (b/2)/ N;

    % ---- Critical Buckling Stress (skin)
    sigma_cr = (k*pi^2*E_skin/(12*(1-nu^2))) ...
                * (t_skin/panel_length)^2;

    % ---- Critical Shear Buckling (web)
    tau_cr = (k*pi^2*E_web/(12*(1-nu^2))) ...
             * (t_web/panel_length)^2;

    % ---- Use MAX stress along span
    sigma_max = max(abs(sigma_bend));
    tau_max   = max(abs(tau_web));

    % ---- Margins
    MS_stringer_buckling(i) = sigma_cr/sigma_max - 1;
    MS_web_shear(i) = tau_cr/tau_max - 1;
    MS_bending_material(i) = sigma_allow/sigma_max - 1;

end
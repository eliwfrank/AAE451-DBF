%% Clark Y 3-cell torsion (symbolic) + carbon tube spars + balsa supports + plots

%% ---------------- USER INPUTS ----------------
csvFile = "ClarkY";    % columns: x_over_c, y_over_c
c = 160.1;                 % [mm]
xCuts = [0.0 0.25 0.75 1.0]*c;    % cell boundaries [mm] => 3 cells

T = T(1)*1000;                  % [N*mm]  (as requested)

% Thicknesses [mm]
t_skin = 0.038;            % Monokote skin thickness [mm] (as requested)
t_web  = 1.50;             % balsa support/web "effective" thickness [mm] (edit if needed)

% Shear moduli [N/mm^2] (EDIT THESE to whatever your team assumes)
% (1 N/mm^2 = 1 MPa)
G_skin = 80;               % Monokote effective shear modulus (very uncertain; placeholder)
G_web  = 180;              % balsa shear modulus (placeholder; depends on density/grain)
G_tube = 5000;            % carbon tube effective torsional shear modulus (placeholder)

% Cell 1 (LE -> front spar) skin = balsa sheeting
t_sheet = 1;      % [mm] 1/16" ~ 1.59 mm (edit to your actual)
G_sheet = 280;      % [N/mm^2] effective balsa shear modulus (edit)

% Carbon tube geometry at x/c = 0.25 and 0.75 (for plot + torsion stiffness)
tubeStations = [0.25 0.75]*c;   % [mm]
tubeOD = 8.0;                   % [mm] outer diameter
tubeID = 6.0;                   % [mm] inner diameter

% Balsa support block for plotting (purely visual)
balsaBlockWidth  = 12;          % [mm]
balsaBlockHeight = 12;          % [mm]
%% --------------------------------------------

%% 1) Load airfoil and scale
M = readmatrix(csvFile);
x = M(:,1)*c;
y = M(:,2)*c;

% Quick sanity check: Clark Y thickness should be ~18-19 mm at c ~160 mm
fprintf("Sanity check: tmax = %.2f mm\n", max(y)-min(y));

% Split into upper/lower by LE (min x)
[~, iLE] = min(x);
x_upper = x(1:iLE);  y_upper = y(1:iLE);
x_lower = x(iLE:end); y_lower = y(iLE:end);

% Ensure both go LE -> TE (increasing x)
if x_upper(1) > x_upper(end), x_upper = flipud(x_upper); y_upper = flipud(y_upper); end
if x_lower(1) > x_lower(end), x_lower = flipud(x_lower); y_lower = flipud(y_lower); end

yu = @(xx) interp1(x_upper, y_upper, xx, "linear", "extrap");
yl = @(xx) interp1(x_lower, y_lower, xx, "linear", "extrap");

% Utility: polyline length
polylen = @(X,Y) sum(sqrt(diff(X).^2 + diff(Y).^2));

% Clamp xCuts to airfoil range
xmin = min(x); xmax = max(x);
xCuts = max(min(xCuts, xmax), xmin);
xCuts = unique(xCuts, "stable");

nCells = numel(xCuts)-1;     % should be 3
assert(nCells==3, "This script expects 3 cells (2 internal webs).");

%% 2) Compute cell areas A1..A3 from polygons
A = zeros(nCells,1);
cells = struct();

for k = 1:nCells
    xa = xCuts(k); xb = xCuts(k+1);

    Xu = x_upper(x_upper>=xa & x_upper<=xb);
    Xu = unique([xa; Xu(:); xb], "stable");
    Yu = yu(Xu);

    Xl = x_lower(x_lower>=xa & x_lower<=xb);
    Xl = unique([xa; Xl(:); xb], "stable");
    Yl = yl(Xl);
    Xl = flipud(Xl); Yl = flipud(Yl);

    % CCW polygon: upper xa->xb, down at xb, lower xb->xa, up at xa
    Xpoly = [Xu(:); xb; Xl(:); xa];
    Ypoly = [Yu(:); yl(xb); Yl(:); yu(xa)];

    A(k) = polyarea(Xpoly, Ypoly);  % mm^2

    cells(k).xa = xa; cells(k).xb = xb;
    cells(k).X = Xpoly; cells(k).Y = Ypoly;
end

fprintf("\nCell areas [mm^2]:\n");
disp(table((1:nCells)', A, 'VariableNames', {'Cell','Area_mm2'}));

%% 3) Compute wall "compliance terms" r = sum(L/(G*t)) for each cell
L_upper = zeros(nCells,1);
L_lower = zeros(nCells,1);

for k = 1:nCells
    xa = xCuts(k); xb = xCuts(k+1);

    Xu = x_upper(x_upper>=xa & x_upper<=xb);
    Xu = unique([xa; Xu(:); xb], "stable");
    Yu = yu(Xu);
    L_upper(k) = polylen(Xu, Yu);

    Xl = x_lower(x_lower>=xa & x_lower<=xb);
    Xl = unique([xa; Xl(:); xb], "stable");
    Yl = yl(Xl);
    L_lower(k) = polylen(Xl, Yl);
end

% Web lengths at internal cuts xCuts(2) and xCuts(3)
xWeb12 = xCuts(2);
xWeb23 = xCuts(3);
L_web12 = abs(yu(xWeb12) - yl(xWeb12));
L_web23 = abs(yu(xWeb23) - yl(xWeb23));

% --- Skin material per cell (Cell 1 = balsa, Cells 2-3 = Monokote) ---
G_skin_cell = [G_sheet; G_skin; G_skin];    % [N/mm^2]
t_skin_cell = [t_sheet; t_skin; t_skin];    % [mm]

% Convert to compliance-like r terms: L/(G*t) per cell
rU = L_upper ./ (G_skin_cell .* t_skin_cell);
rL = L_lower ./ (G_skin_cell .* t_skin_cell);

% Webs (still balsa/ply supports)
rW12 = L_web12/(G_web*t_web);
rW23 = L_web23/(G_web*t_web);

disp(table((1:3)', G_skin_cell, t_skin_cell, ...
    'VariableNames', {'Cell','G_skin_used','t_skin_used'}));

fprintf("\nWall summary:\n");
disp(table((1:nCells)', L_upper, L_lower, rU, rL, 'VariableNames', ...
    {'Cell','L_upper_mm','L_lower_mm','r_upper=L/(G*t)','r_lower=L/(G*t)'}));
disp(table(["web@0.25c";"web@0.75c"], [L_web12; L_web23], [rW12; rW23], ...
    'VariableNames', {'Web','L_web_mm','r_web=L/(G*t)'}));

%% 4) Carbon tube torsional stiffness sum(GJ) added to torque equation
Ro = tubeOD/2; Ri = tubeID/2;
J_tube = (pi/2)*(Ro^4 - Ri^4);         % mm^4
GJ_each = G_tube * J_tube;             % N*mm^2
sumGJ = numel(tubeStations) * GJ_each;

fprintf("\nTube: J = %.3e mm^4, each GJ = %.3e N*mm^2, sum(GJ) = %.3e N*mm^2\n", ...
    J_tube, GJ_each, sumGJ);

%% 5) SYMBOLIC solve (professor style) for 3 cells
syms q1 q2 q3 theta1 theta2 theta3

A1 = A(1); A2 = A(2); A3 = A(3);

% Compatibility for each cell:
% sum_over_walls( q_wall * L/(G*t) ) - 2*A_i*theta_i = 0
eq1 = q1*(rU(1)+rL(1)) + (q1-q2)*rW12 - 2*A1*theta1 == 0;
eq2 = q2*(rU(2)+rL(2)) + (q2-q1)*rW12 + (q2-q3)*rW23 - 2*A2*theta2 == 0;
eq3 = q3*(rU(3)+rL(3)) + (q3-q2)*rW23 - 2*A3*theta3 == 0;

% Torque equilibrium includes tube torque: 2*sum(q_i*A_i) + sum(GJ)*theta = T
eq4 = 2*q1*A1 + 2*q2*A2 + 2*q3*A3 + sumGJ*theta1 == T;

% Same twist rate for all cells
eq5 = theta1 - theta2 == 0;
eq6 = theta2 - theta3 == 0;

sol = solve([eq1 eq2 eq3 eq4 eq5 eq6], [q1 q2 q3 theta1 theta2 theta3]);

Q1 = double(vpa(sol.q1,10));
Q2 = double(vpa(sol.q2,10));
Q3 = double(vpa(sol.q3,10));
th = double(vpa(sol.theta1,12));

fprintf("\n=== SYMBOLIC SOLUTION ===\n");
fprintf("q1 = %.6f N/mm\n", Q1);
fprintf("q2 = %.6f N/mm\n", Q2);
fprintf("q3 = %.6f N/mm\n", Q3);
fprintf("theta' = %.8e rad/mm\n", th);

%% 6) Plot airfoil + balsa webs + carbon tube circles (geometry view)
figure; hold on; grid on; axis equal;

plot(x/c, y/c, 'k-', 'LineWidth', 1.5);

% Web lines (balsa supports) at x/c = 0.25 and 0.75
webXs = [xWeb12 xWeb23];
for k = 1:numel(webXs)
    xx = webXs(k);
    plot([xx xx]/c, [yl(xx) yu(xx)]/c, 'b-', 'LineWidth', 2);
end

% Balsa blocks around tube centers (visual only)
for k = 1:numel(tubeStations)
    xc = tubeStations(k);
    yc = 0.5*(yu(xc)+yl(xc));  % mid-thickness

    x0 = xc - balsaBlockWidth/2;  x1r = xc + balsaBlockWidth/2;
    y0 = yc - balsaBlockHeight/2; y1r = yc + balsaBlockHeight/2;

    plot([x0 x1r x1r x0 x0]/c, [y0 y0 y1r y1r y0]/c, 'Color', [0.1 0.6 0.1], 'LineWidth', 1.5);
end

% Carbon tube circles
thang = linspace(0,2*pi,200);
for k = 1:numel(tubeStations)
    xc = tubeStations(k);
    yc = 0.5*(yu(xc)+yl(xc));  % mid-thickness point (adjust if you want)

    rO = tubeOD/2; rI = tubeID/2;
    plot((xc + rO*cos(thang))/c, (yc + rO*sin(thang))/c, 'r-', 'LineWidth', 1.5);
    plot((xc + rI*cos(thang))/c, (yc + rI*sin(thang))/c, 'r--', 'LineWidth', 1.0);
end

% Label cells
for k = 1:nCells
    cx = mean(cells(k).X)/c;
    cy = mean(cells(k).Y)/c;
    text(cx, cy, sprintf("Cell %d", k), 'FontWeight','bold');
end

xlabel('x/c'); ylabel('y/c');
title(sprintf('Clark Y (c=%.1f mm), 3 cells @ 0.25c & 0.75c, T=%.0f Nmm', c, T));
legend({'Airfoil','Balsa webs','Balsa blocks','Tube OD','Tube ID'}, 'Location','bestoutside');

%% 7) Shear-flow colored plot that FOLLOWS the airfoil lines
% Wall shear flows:
q_upper = [Q1 Q2 Q3];     % upper skin in each cell
q_lower = [Q1 Q2 Q3];     % lower skin in each cell
q_web12 = Q1 - Q2;        % web between cell 1 and 2
q_web23 = Q2 - Q3;        % web between cell 2 and 3

% For color scaling, collect magnitudes
q_all = [abs(q_upper) abs(q_lower) abs(q_web12) abs(q_web23)];
qmin = min(q_all);
qmax = max(q_all);

cmap = parula(256);

mapIdx = @(val) max(1, min(256, round(1 + 255*(val-qmin)/(qmax-qmin+eps))));

figure; hold on; grid on; axis equal;
plot(x/c, y/c, 'k-', 'LineWidth', 1.2);  % base airfoil in black

% Helper to get polyline segment points on upper/lower between xa and xb
getSeg = @(xSurf, ySurf, xa, xb, f) deal( ...
    unique([xa; xSurf(xSurf>=xa & xSurf<=xb); xb], "stable"), ...
    f(unique([xa; xSurf(xSurf>=xa & xSurf<=xb); xb], "stable")) );

% Plot colored skins (upper + lower) for each cell bay
for k = 1:nCells
    xa = xCuts(k); xb = xCuts(k+1);

    % Upper surface segment points
    [Xu, Yu] = getSeg(x_upper, y_upper, xa, xb, yu);
    colU = cmap(mapIdx(abs(q_upper(k))), :);
    plot(Xu/c, Yu/c, '-', 'LineWidth', 5, 'Color', colU);

    % Lower surface segment points
    [Xl, Yl] = getSeg(x_lower, y_lower, xa, xb, yl);
    colL = cmap(mapIdx(abs(q_lower(k))), :);
    plot(Xl/c, Yl/c, '-', 'LineWidth', 5, 'Color', colL);

    % Label near middle of each segment
    midU = round(numel(Xu)/2);
    text(Xu(midU)/c, Yu(midU)/c, sprintf('q%d=%.2f', k, q_upper(k)), ...
        'FontSize', 8, 'BackgroundColor', 'w', 'Margin', 1);

    midL = round(numel(Xl)/2);
    text(Xl(midL)/c, Yl(midL)/c, sprintf('q%d=%.2f', k, q_lower(k)), ...
        'FontSize', 8, 'BackgroundColor', 'w', 'Margin', 1);
end

% Plot colored webs
% Web at 0.25c
colW12 = cmap(mapIdx(abs(q_web12)), :);
plot([xWeb12 xWeb12]/c, [yl(xWeb12) yu(xWeb12)]/c, '-', 'LineWidth', 6, 'Color', colW12);
text(xWeb12/c, 0.5*(yl(xWeb12)+yu(xWeb12))/c, sprintf('q12=%.2f', q_web12), ...
    'FontSize', 8, 'BackgroundColor', 'w', 'Margin', 1);

% Web at 0.75c
colW23 = cmap(mapIdx(abs(q_web23)), :);
plot([xWeb23 xWeb23]/c, [yl(xWeb23) yu(xWeb23)]/c, '-', 'LineWidth', 6, 'Color', colW23);
text(xWeb23/c, 0.5*(yl(xWeb23)+yu(xWeb23))/c, sprintf('q23=%.2f', q_web23), ...
    'FontSize', 8, 'BackgroundColor', 'w', 'Margin', 1);

% Add tubes and balsa blocks lightly for context (optional)
for k = 1:numel(tubeStations)
    xc = tubeStations(k);
    yc = 0.5*(yu(xc)+yl(xc));
    rO = tubeOD/2; rI = tubeID/2;
    plot((xc + rO*cos(thang))/c, (yc + rO*sin(thang))/c, 'k-', 'LineWidth', 1.0);
    plot((xc + rI*cos(thang))/c, (yc + rI*sin(thang))/c, 'k--', 'LineWidth', 0.8);
end

% --- Fix colorbar scaling to match q-values ---
colormap(cmap);          % use the same colormap you used for line colors
clim([qmin qmax]);      % force color axis to the actual |q| range
cb = colorbar;
ylabel(cb, '|q| [N/mm]');
% ---------------------------------------------

xlabel('x/c'); ylabel('y/c');
title(sprintf('Torsional Shear Flow (T = %.1f N-mm) - With D-Box', T));


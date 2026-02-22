%% Clark Y 2-cell torsion (symbolic) + carbon tube spars + balsa supports + plots
% Cells:
%   Cell 1: LE D-box        from 0    to 0.25c
%   Cell 2: Main wing box   from 0.25c to 0.75c
% Middle cell is CLOSED by a rear web at x/c = 0.75 (IMPORTANT).
% NO trailing-edge (TE) cell is modeled.

%% ---------------- USER INPUTS ----------------
csvFile = "NACA0012.csv";      % columns: x_over_c, y_over_c
chord = c_tail_h*1000;                   % [mm]

T_ht = T_h(1)*1000; % in N*mm

% Cell boundaries (mm) => 2 cells (NO TE cell)
xCuts = [0.0 0.25 0.6]*chord;   % [mm]

% Thicknesses [mm]
t_skin = 0.038;              % Monokote skin thickness [mm]
t_web  = 1.00;               % balsa support/web "effective" thickness [mm]

% Shear moduli [N/mm^2] (EDIT as needed)
G_skin  = 80;                % Monokote effective shear modulus (placeholder)
G_web   = 180;               % balsa shear modulus (placeholder)
G_tube  = 5000;              % carbon tube effective shear modulus (placeholder)

% Cell 1 (LE -> front spar) skin = balsa sheeting (D-box)
t_sheet = 1.00;              % [mm] balsa sheet thickness
G_sheet = 180;               % [N/mm^2] balsa shear modulus (placeholder)

% --- Spar types (front stays carbon tube, rear becomes balsa I-beam) ---
sparStations = [0.25 0.6]*chord;     % [mm] (same as your tubeStations)
tubeStations = sparStations;

% Carbon tube (front spar only)
tubeOD = 5;                 % [mm]
tubeID = 4;                 % [mm]

% Balsa I-beam (rear spar at x/c = 0.6)
b_fl_I   = 12.0;             % [mm] flange width
t_fl_I   = 1.5;              % [mm] flange thickness
t_web_I  = 2.0;              % [mm] web thickness
fill_I   = 0.90;             % [-] fraction of local airfoil thickness used as I-beam depth
G_Ibeam  = G_web;            % [N/mm^2] use balsa shear modulus (or define separately)

% Balsa support block for plotting (visual only)
balsaBlockWidth  = 12;        % [mm]
balsaBlockHeight = 12;        % [mm]
%% --------------------------------------------

%% 1) Load airfoil and scale
M = readmatrix(csvFile);
x = M(:,1)*chord;
y = M(:,2)*chord;

% Split into upper/lower by LE (min x)
[~, iLE] = min(x);
x_upper = x(1:iLE);   y_upper = y(1:iLE);
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

nCells = numel(xCuts)-1;     % should be 2
assert(nCells==2, "This script expects 2 cells (D-box + main box).");

%% 2) Compute cell areas A1..A2 from polygons
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

    % CCW polygon
    Xpoly = [Xu(:); xb; Xl(:); xa];
    Ypoly = [Yu(:); yl(xb); Yl(:); yu(xa)];

    A(k) = polyarea(Xpoly, Ypoly);  % mm^2
    cells(k).xa = xa; cells(k).xb = xb;
    cells(k).X = Xpoly; cells(k).Y = Ypoly;
end

fprintf("\nCell areas [mm^2]:\n");
disp(table((1:nCells)', A, 'VariableNames', {'Cell','Area_mm2'}));

%% 3) Compute wall compliance terms r = sum(L/(G*t)) for each cell
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

% Internal web at x/c = 0.25 (between cell 1 and cell 2)
xWeb12  = xCuts(2);
L_web12 = abs(yu(xWeb12) - yl(xWeb12));
rW12    = L_web12/(G_web*t_web);

% Rear closing web at x/c = 0.75 (closes cell 2)
xWeb2   = xCuts(3);
L_web2  = abs(yu(xWeb2) - yl(xWeb2));
rW2     = L_web2/(G_web*t_web);

% Skin material per cell (Cell 1 = balsa sheet, Cell 2 = Monokote)
G_skin_cell = [G_sheet; G_skin];    % [N/mm^2]
t_skin_cell = [t_sheet; t_skin];    % [mm]

% Compliance-like r terms
rU = L_upper ./ (G_skin_cell .* t_skin_cell);
rL = L_lower ./ (G_skin_cell .* t_skin_cell);

fprintf("\nWall summary:\n");
disp(table((1:nCells)', L_upper, L_lower, rU, rL, ...
    'VariableNames', {'Cell','L_upper_mm','L_lower_mm','rU=L/(G*t)','rL=L/(G*t)'}));
disp(table(["web@0.25c";"web@0.75c"], [L_web12; L_web2], [rW12; rW2], ...
    'VariableNames', {'Web','L_web_mm','rW=L/(G*t)'}));

%% 4) Spar torsional stiffness sum(GJ) added to torque equation
% Front spar = carbon tube @ x/c=0.25
Ro = tubeOD/2; Ri = tubeID/2;
J_tube = (pi/2)*(Ro^4 - Ri^4);     % [mm^4]
GJ_tube = G_tube * J_tube;         % [N*mm^2]

% Rear spar = balsa I-beam @ x/c=0.6 (open-section torsion approx)
xI = tubeStations(2);                              % [mm] station at 0.6c
h_local = abs(yu(xI) - yl(xI));                    % [mm] local airfoil thickness
h_beam  = fill_I * h_local;                        % [mm] chosen I-beam depth
h_web   = max(h_beam - 2*t_fl_I, 0);               % [mm] clear web height

% St-Venant torsion constant for open thin-walled sections (engineering approx):
% J ≈ (1/3) * Σ (b_i * t_i^3)
J_Ibeam = (1/3) * (2*b_fl_I*(t_fl_I^3) + h_web*(t_web_I^3));   % [mm^4]
GJ_Ibeam = G_Ibeam * J_Ibeam;                                  % [N*mm^2]

% Total added stiffness
sumGJ = GJ_tube + GJ_Ibeam;

fprintf("\nSpars torsion stiffness:\n");
fprintf("Tube:  J = %.3e mm^4, GJ = %.3e N*mm^2\n", J_tube,  GJ_tube);
fprintf("Ibeam: J = %.3e mm^4, GJ = %.3e N*mm^2 (h_beam=%.2f mm)\n", J_Ibeam, GJ_Ibeam, h_beam);
fprintf("sum(GJ) = %.3e N*mm^2\n", sumGJ);

%% 5) SYMBOLIC solve (2 closed cells: D-box + main box)
syms q1 q2 theta1 theta2
A1 = A(1); A2 = A(2);

% Compatibility:
% Cell 1: upper+lower skins + shared web (0.25c)
eq1 = q1*(rU(1)+rL(1)) + (q1-q2)*rW12 - 2*A1*theta1 == 0;

% Cell 2: upper+lower skins + shared web (0.25c) + rear closing web (0.75c)
eq2 = q2*(rU(2)+rL(2)) + (q2-q1)*rW12 + q2*rW2 - 2*A2*theta2 == 0;

% Torque equilibrium includes tube torque
eq3 = 2*q1*A1 + 2*q2*A2 + sumGJ*theta1 == T_ht;

% Same twist rate
eq4 = theta1 - theta2 == 0;

sol = solve([eq1 eq2 eq3 eq4], [q1 q2 theta1 theta2]);

Q1 = double(vpa(sol.q1,10));
Q2 = double(vpa(sol.q2,10));
th = double(vpa(sol.theta1,12));

fprintf("\n=== SYMBOLIC SOLUTION (2 CLOSED cells) ===\n");
fprintf("q1 = %.6f N/mm\n", Q1);
fprintf("q2 = %.6f N/mm\n", Q2);
fprintf("theta' = %.8e rad/mm\n", th);

%% 6) Geometry view: airfoil + webs + balsa blocks + carbon tubes

thang = linspace(0,2*pi,200);

%% 7) Shear-flow colored plot following airfoil lines (2 cells only)
% Wall shear flows:
q_upper = [Q1 Q2];        % upper skin in each cell bay
q_lower = [Q1 Q2];        % lower skin in each cell bay
q_web12 = Q1 - Q2;        % shared web @0.25c
q_web2  = Q2;             % rear closing web @0.75c (outside is open)

q_all = [abs(q_upper) abs(q_lower) abs(q_web12) abs(q_web2)];
qmin = min(q_all);
qmax = max(q_all);

cmap = parula(256);
mapIdx = @(val) max(1, min(256, round(1 + 255*(val-qmin)/(qmax-qmin+eps))));

figure; hold on; grid on; axis equal;
plot(x/chord, y/chord, 'k-', 'LineWidth', 1.2);

% segment helper
getSeg = @(xSurf, xa, xb, f) deal( ...
    unique([xa; xSurf(xSurf>=xa & xSurf<=xb); xb], "stable"), ...
    f(unique([xa; xSurf(xSurf>=xa & xSurf<=xb); xb], "stable")) );

for k = 1:nCells
    xa = xCuts(k); xb = xCuts(k+1);

    [Xu, Yu] = getSeg(x_upper, xa, xb, yu);
    colU = cmap(mapIdx(abs(q_upper(k))), :);
    plot(Xu/chord, Yu/chord, '-', 'LineWidth', 5, 'Color', colU);

    [Xl, Yl] = getSeg(x_lower, xa, xb, yl);
    colL = cmap(mapIdx(abs(q_lower(k))), :);
    plot(Xl/chord, Yl/chord, '-', 'LineWidth', 5, 'Color', colL);

    midU = round(numel(Xu)/2);
    text(Xu(midU)/chord, Yu(midU)/chord, sprintf('q%d=%.2f', k, q_upper(k)), ...
        'FontSize', 8, 'BackgroundColor', 'w', 'Margin', 1);

    midL = round(numel(Xl)/2);
    text(Xl(midL)/chord, Yl(midL)/chord, sprintf('q%d=%.2f', k, q_lower(k)), ...
        'FontSize', 8, 'BackgroundColor', 'w', 'Margin', 1);
end

% web @ 0.25c
colW12 = cmap(mapIdx(abs(q_web12)), :);
plot([xWeb12 xWeb12]/chord, [yl(xWeb12) yu(xWeb12)]/chord, '-', 'LineWidth', 6, 'Color', colW12);
text(xWeb12/chord, 0.5*(yl(xWeb12)+yu(xWeb12))/chord, sprintf('q12=%.2f', q_web12), ...
    'FontSize', 8, 'BackgroundColor', 'w', 'Margin', 1);

% web @ 0.75c (rear closing web)
colW2 = cmap(mapIdx(abs(q_web2)), :);
plot([xWeb2 xWeb2]/chord, [yl(xWeb2) yu(xWeb2)]/chord, '-', 'LineWidth', 6, 'Color', colW2);
text(xWeb2/chord, 0.5*(yl(xWeb2)+yu(xWeb2))/chord, sprintf('q2=%.2f', q_web2), ...
    'FontSize', 8, 'BackgroundColor', 'w', 'Margin', 1);

% tubes for context

xc = tubeStations(1);
yc = 0.5*(yu(xc)+yl(xc));
rO = tubeOD/2; rI = tubeID/2;
plot((xc + rO*cos(thang))/chord, (yc + rO*sin(thang))/chord, 'k-', 'LineWidth', 1.0);
plot((xc + rI*cos(thang))/chord, (yc + rI*sin(thang))/chord, 'k--', 'LineWidth', 0.8);


colormap(cmap);
clim([qmin qmax]);
cb = colorbar;
ylabel(cb, '|q| [N/mm]');

xlabel('x/c'); ylabel('y/c');
title(sprintf('Horz. Stab. Torsional Shear Flow (T = %.1f N-mm)', T_ht));

%% Flexural shear flow (MIDDLE CELL ONLY) + q0 + Shear Center (x_sc)
% Clark-Y, spars at x/c = 0.25 and 0.75
% Middle cell only: xF -> xR with 3 stringers top and 2 bottom
%
% Computes:
%   - booms (spars + stringers in middle bay only)
%   - Iy from boom idealization
%   - basic shear flow q'(s) via boom-walk
%   - q0 from zero-twist compatibility: integral(q/(tG) ds)=0
%   - shear center x_sc from M = ∮ q (x dz - z dx),  x_sc = M/Vz
%   - plot q_total = q' + q0 on the middle-cell perimeter (x/c, z/c)
%
% Notes:
% - This is an engineering approximation for a single closed cell.
% - Uses skin/web thickness + shear modulus only for q0 (compatibility), not for q'.


%% ========== USER INPUTS ==========
csvFile = "NACA0012.csv";     % airfoil file: columns [x/c, z/c]
chord = c_tail_h*1000;                  % chord [mm]

Vz_ht = V_t(1);

xF = 0.25;                  % front spar x/c
xR = 0.6;                  % rear  spar x/c

% ----- Materials for q0 (compatibility) -----
% Use your best estimates (placeholders OK for now)
t_skin = 0.038;             % skin thickness [mm] (Monokote)
G_skin = 80;                % [N/mm^2] effective shear modulus (Monokote placeholder)

t_web  = 1.0;               % web thickness [mm]
G_web  = 180;               % [N/mm^2] web shear modulus (balsa/ply placeholder)

% ----- Carbon tube spar boom area -----
Do = 5;                     % [mm] OD
Di = 4;                     % [mm] ID
A_tube = (pi/4)*(Do^2 - Di^2);  % [mm^2]

% ----- Rear spar (at xR) as BALSA I-beam (boom-equivalent) -----
b_fl_R   = 10.0;    % [mm] flange width
t_fl_R   = 1.5;     % [mm] flange thickness
t_web_R  = 2.0;     % [mm] web thickness
fill_R   = 0.90;    % [-] fraction of local thickness used as I-beam depth (keep < 1)
webSplit = 0.50;    % [-] fraction of web area to lump into each flange-boom (0.5 is common quick approx)

% ----- Balsa stringers (boom areas) -----
b_str = 3.0;                % [mm] stringer width
k_str = 0.30;               % height = 30% of local thickness
hmax_str = 6.0;             % clamp max height [mm]

% Middle cell stringers only (3 top, 2 bottom)
x_str_top = [0.32 0.42 0.52];
x_str_bot = [0.37 0.47];

% Plot controls
showLabels = true;
Nsamp_plot = 120;           % smooth plotting of each panel along airfoil
Nsamp_int  = 200;           % smooth integration for q0 + shear center

%% ========== 1) READ AIRFOIL / BUILD INTERPOLANTS ==========
M = readmatrix(csvFile);
x_all = M(:,1); z_all = M(:,2);
mask = ~(isnan(x_all) | isnan(z_all));
x_all = x_all(mask); z_all = z_all(mask);

% Split into upper/lower by LE index (min x)
[~, iLE] = min(x_all);
x1 = x_all(1:iLE);   z1 = z_all(1:iLE);
x2 = x_all(iLE:end); z2 = z_all(iLE:end);

% Decide which is upper by mean z
if mean(z1) > mean(z2)
    xU = x1; zU = z1;  xL = x2; zL = z2;
else
    xU = x2; zU = z2;  xL = x1; zL = z1;
end

% Unique + sort increasing x for interpolation
[xU, iu] = unique(xU, 'stable'); zU = zU(iu);
[xL, il] = unique(xL, 'stable'); zL = zL(il);
[xU, su] = sort(xU); zU = zU(su);
[xL, sl] = sort(xL); zL = zL(sl);

zu_c = @(xc) interp1(xU, zU, xc, 'linear', 'extrap');  % z/c on upper
zl_c = @(xc) interp1(xL, zL, xc, 'linear', 'extrap');  % z/c on lower

%% ========== 2) BUILD BOOMS (MIDDLE CELL ONLY) ==========
% Front Spar
t_avail = @(xc) (zu_c(xc) - zl_c(xc))*chord;  % local thickness [mm]
A_balsa = @(xc) b_str * min(k_str*t_avail(xc), hmax_str);

% ---- Rear spar I-beam: compute equivalent boom areas at xR ----
h_beam_R = fill_R * t_avail(xR);                 % [mm] overall I-beam depth (fits inside section)
h_web_R  = max(h_beam_R - 2*t_fl_R, 0);          % [mm] clear web height

A_fl_R   = b_fl_R * t_fl_R;                      % [mm^2] one flange area
A_web_R  = t_web_R * h_web_R;                    % [mm^2] web area

% Boom-equivalent areas (dominant bending carried by flanges; optionally lump some web area)
B_R_top  = A_fl_R + webSplit*A_web_R;            % [mm^2]
B_R_bot  = A_fl_R + webSplit*A_web_R;            % [mm^2]

% (Optional) sanity print
fprintf("Rear spar I-beam: h=%.2f mm, A_fl=%.2f mm^2, A_web=%.2f mm^2, Btop=Bbot=%.2f mm^2\n", ...
        h_beam_R, A_fl_R, A_web_R, B_R_top);

booms = struct('name',{},'xc',{},'x_mm',{},'z_mm',{},'B',{});
addBoom = @(name,xc,zmm,B) struct('name',string(name),'xc',xc,'x_mm',xc*chord,'z_mm',zmm,'B',B);

% Spar booms at xF and xR (top/bot)
booms(end+1) = addBoom("F_spar_top", xF, zu_c(xF)*chord, A_tube);
booms(end+1) = addBoom("F_spar_bot", xF, zl_c(xF)*chord, A_tube);
booms(end+1) = addBoom("R_spar_top", xR, zu_c(xR)*chord, B_R_top);
booms(end+1) = addBoom("R_spar_bot", xR, zl_c(xR)*chord, B_R_bot);

% Middle-cell stringers only
for k = 1:numel(x_str_top)
    xc = x_str_top(k);
    booms(end+1) = addBoom("StrTop_"+k, xc, zu_c(xc)*chord, A_balsa(xc));
end
for k = 1:numel(x_str_bot)
    xc = x_str_bot(k);
    booms(end+1) = addBoom("StrBot_"+k, xc, zl_c(xc)*chord, A_balsa(xc));
end

%% ========== 3) Iy (boom idealization) ==========
Bvec = [booms.B]';
zvec = [booms.z_mm]';
zbar = sum(Bvec .* zvec) / sum(Bvec);     % centroid in z
zc   = zvec - zbar;
Iy   = sum(Bvec .* (zc.^2));              % [mm^4]
fprintf("Iy (boom idealization, middle cell) = %.6e mm^4\n", Iy);

%% ========== 4) BUILD STATIONS AROUND MIDDLE CELL ==========
stations = buildStationsCell_MID_ONLY(booms, xF, xR, zu_c, zl_c, chord, zbar);

%% ========== 5) BASIC q'(s) BY BOOM-WALK ==========
[qPanels, ~] = shearFlowPanelsFromStations(stations, Vz_ht, Iy);

%% ========== 6) SOLVE q0 (ZERO TWIST) + SHEAR CENTER x_sc ==========
[q0, x_sc_mm, x_sc_over_c, theta_check] = solve_q0_and_shearcenter( ...
    stations, qPanels, Vz_ht, xF, xR, zu_c, zl_c, chord, Nsamp_int, t_skin, G_skin, t_web, G_web);

fprintf("\n--- Closed-cell correction (middle cell) ---\n");
fprintf("q0 = %.6f N/mm\n", q0);
fprintf("x_sc = %.3f mm  (x_sc/c = %.4f)\n", x_sc_mm, x_sc_over_c);
fprintf("theta' check (should be ~0): %.3e rad/mm\n", theta_check);

% Total shear flow on each panel
qPanels_tot = qPanels + q0;

%% ========== 7) PLOT TOTAL q = q' + q0 ==========
qAbs = abs(qPanels_tot(:));
qMin = min(qAbs); qMax = max(qAbs);
if abs(qMax-qMin) < 1e-12, qMax = qMin + 1; end

figure; hold on; grid on; axis equal;
title(sprintf("Horiz. Stab. Flexural Shear Flow with q0 (Vz = %.1f N)", Vz_ht));
xlabel("x/c"); ylabel("z/c");

% Airfoil outline
xc_plot = linspace(0,1,600);
plot(xc_plot, zu_c(xc_plot), 'k-', 'LineWidth', 1.2);
plot(xc_plot, zl_c(xc_plot), 'k-', 'LineWidth', 1.2);

% Webs (front + rear)
plot([xF xF], [zu_c(xF) zl_c(xF)], 'k-', 'LineWidth', 2.0);
plot([xR xR], [zu_c(xR) zl_c(xR)], 'k-', 'LineWidth', 2.0);

% Booms
scatter([booms.xc], [booms.z_mm]/chord, 35, 'filled');
for i=1:numel(booms)
    text(booms(i).xc, booms(i).z_mm/chord, "  "+booms(i).name, 'FontSize', 8);
end

% Mark shear center on x-axis (z=0 line) for reference
plot([x_sc_over_c x_sc_over_c], [-0.02 0.02], 'r-', 'LineWidth', 2);
text(x_sc_over_c, 0.025, sprintf("x_sc/c=%.3f", x_sc_over_c), 'Color','r', 'FontWeight','bold');

colormap(parula(256));
set(gca,'CLim',[qMin qMax]);

% Dummy mappable for colorbar stability
hDummy = scatter(nan,nan,1,nan,'filled'); 
set(hDummy,'Visible','off');

% Colored panels (total)
plotColoredPanelsAirfoil_xc(stations, qPanels_tot, qMin, qMax, showLabels, chord, zu_c, zl_c, Nsamp_plot);

cb = colorbar;
cb.Label.String = "|q| (N/mm)";

T_h = V_t*abs(c_tail_h/4 - x_sc_over_c*c_tail_h); % will be in N*m, this is Wing Torsion

%% ===================== LOCAL FUNCTIONS =====================

function stations = buildStationsCell_MID_ONLY(booms, xF, xR, zu_c, zl_c, c, zbar)
% Middle cell perimeter:
% F_top -> upper to R_top -> rear web down -> lower back -> front web up

tol = 1e-9;
top = booms([booms.z_mm] > 0 & [booms.xc] >= xF-tol & [booms.xc] <= xR+tol);
bot = booms([booms.z_mm] < 0 & [booms.xc] >= xF-tol & [booms.xc] <= xR+tol);

[~,idx] = sort([top.xc],'ascend');  top = top(idx);     % xF -> xR
[~,idx] = sort([bot.xc],'descend'); bot = bot(idx);     % xR -> xF

stations = struct('name',{},'x_mm',{},'z_mm',{},'B',{},'zc',{},'isSkin',{},'side',{});

for i=1:numel(top)
    stations(end+1) = makeStationFromBoom(top(i), "top", zbar);
end

% rear web down (to lower at xR)
stations(end+1) = struct('name',"R_web_bot", 'x_mm',xR*c, 'z_mm',zl_c(xR)*c, 'B',0, 'zc',0, ...
                         'isSkin', false, 'side',"web");

for i=1:numel(bot)
    stations(end+1) = makeStationFromBoom(bot(i), "bot", zbar);
end

% front web up closure
stations(end+1) = struct('name',"F_web_top", 'x_mm',xF*c, 'z_mm',zu_c(xF)*c, 'B',0, 'zc',0, ...
                         'isSkin', false, 'side',"web");
end

function st = makeStationFromBoom(boom, side, zbar)
st = struct('name', boom.name, ...
            'x_mm', boom.x_mm, ...
            'z_mm', boom.z_mm, ...
            'B', boom.B, ...
            'zc', boom.z_mm - zbar, ...
            'isSkin', true, ...
            'side', side);
end

function [qPanels, panels] = shearFlowPanelsFromStations(stations, Vz, Iy)
n = numel(stations);
qPanels = zeros(n-1,1);
panels  = zeros(n-1,2,2);

q = 0;  % arbitrary cut reference
for k = 1:(n-1)
    q = q + (Vz/Iy) * (stations(k).B * stations(k).zc);
    x1 = stations(k).x_mm;   z1 = stations(k).z_mm;
    x2 = stations(k+1).x_mm; z2 = stations(k+1).z_mm;
    panels(k,:,:) = [x1 z1; x2 z2];
    qPanels(k) = q;
end
end

function plotColoredPanelsAirfoil_xc(stations, qPanels, qMin, qMax, showLabels, c_mm, zu_c, zl_c, Nsamp)
for k = 1:(numel(stations)-1)
    qv = abs(qPanels(k));

    s1 = stations(k);
    s2 = stations(k+1);

    x1c = s1.x_mm / c_mm;  z1c = s1.z_mm / c_mm;
    x2c = s2.x_mm / c_mm;  z2c = s2.z_mm / c_mm;

    if hypot(x2c-x1c, z2c-z1c) < 1e-10
        continue;
    end

    if s1.isSkin && s2.isSkin && (s1.side=="top" || s1.side=="bot")
        xcs = linspace(x1c, x2c, Nsamp);
        if s1.side=="top"
            zcs = zu_c(xcs);
        else
            zcs = zl_c(xcs);
        end
    else
        xcs = [x1c x2c];
        zcs = [z1c z2c];
    end

    % surface-edge colored line (no polygon closure)
    X = [xcs; xcs];
    Y = [zcs; zcs];
    Z = zeros(size(X));
    C = qv * ones(size(X));
    h = surface(X, Y, Z, C, 'FaceColor','none', 'EdgeColor','interp', 'LineWidth',4);
    set(h,'CDataMapping','scaled');

    if showLabels
        text(mean(xcs), mean(zcs), sprintf(" %.2f", qPanels(k)), 'FontSize', 8);
    end
end
set(gca,'CLim',[qMin qMax]);
end

function [q0, x_sc_mm, x_sc_over_c, theta_check] = solve_q0_and_shearcenter( ...
    stations, qPanels, Vz, xF, xR, zu_c, zl_c, c, Nsamp, t_skin, G_skin, t_web, G_web)
% q0 from zero twist: ∮(q' + q0)/(tG) ds = 0
% shear center (about LE): x_sc = [∮ q(x dz - z dx)] / Vz

    % Cell area A (mm^2) between xF and xR
    xc = linspace(xF, xR, 1200);
    Xu = xc*c;  Zu = zu_c(xc)*c;
    Xl = fliplr(Xu);  Zl = fliplr(zl_c(xc)*c);

    Xpoly = [Xu, Xu(end), Xl, Xu(1)];
    Zpoly = [Zu, Zl(1),  Zl, Zu(1)];
    A_cell = polyarea(Xpoly, Zpoly);  % mm^2

    I1 = 0; I0 = 0;  % for q0
    nSeg = numel(stations)-1;

    % First pass: build I1, I0
    for k = 1:nSeg
        s1 = stations(k);
        s2 = stations(k+1);

        [Xmm, Zmm, isWeb] = panelPolyline_mm(s1, s2, zu_c, zl_c, c, Nsamp);
        if numel(Xmm) < 2, continue; end

        ds = sqrt(diff(Xmm).^2 + diff(Zmm).^2);
        Lk = sum(ds);

        if isWeb
            tk = t_web;  Gk = G_web;
        else
            tk = t_skin; Gk = G_skin;
        end

        I1 = I1 + qPanels(k) * (Lk/(tk*Gk));
        I0 = I0 + (Lk/(tk*Gk));
    end

    q0 = -I1 / I0;
    theta_check = (1/(2*A_cell)) * (I1 + q0*I0);  % should be ~0

    % Second pass: torque integral for shear center
    M = 0;  % N*mm
    for k = 1:nSeg
        s1 = stations(k);
        s2 = stations(k+1);

        [Xmm, Zmm, ~] = panelPolyline_mm(s1, s2, zu_c, zl_c, c, Nsamp);
        if numel(Xmm) < 2, continue; end

        % J = ∮(x dz - z dx) along the panel polyline (piecewise linear)
        J = 0;
        for j = 1:(numel(Xmm)-1)
            x1 = Xmm(j);   z1 = Zmm(j);
            x2 = Xmm(j+1); z2 = Zmm(j+1);
            J = J + (x1*z2 - z1*x2);   % mm^2
        end

        qk = qPanels(k) + q0;          % N/mm
        M  = M + qk * J;               % N*mm
    end

    x_sc_mm = -M / Vz;
    x_sc_over_c = x_sc_mm / c;
end

function [Xmm, Zmm, isWeb] = panelPolyline_mm(s1, s2, zu_c, zl_c, c, Nsamp)
% Polyline points for a panel:
% - skin panels follow airfoil (top or bottom)
% - web panels are straight

    x1c = s1.x_mm / c;
    z1c = s1.z_mm / c;
    x2c = s2.x_mm / c;
    z2c = s2.z_mm / c;

    if hypot(x2c-x1c, z2c-z1c) < 1e-12
        Xmm = []; Zmm = []; isWeb = true;
        return;
    end

    isWeb = ~(s1.isSkin && s2.isSkin);

    if ~isWeb
        xcs = linspace(x1c, x2c, Nsamp);
        if string(s1.side) == "top"
            zcs = zu_c(xcs);
        else
            zcs = zl_c(xcs);
        end
        Xmm = xcs * c;
        Zmm = zcs * c;
    else
        Xmm = [x1c x2c] * c;
        Zmm = [z1c z2c] * c;
    end
end
%% Idealized Wingbox Flexural Shear Flow (Boom/Web Method) - FIXED PLOTTING
% Clark-Y airfoil, chord known, spars at x/c=0.25 and 0.75
% Computes:
%  - boom areas (carbon tube + balsa stringers sized from local thickness)
%  - Iy (boom idealization)
%  - basic flexural shear flow q' around each cell (LE cell + middle cell)
%  - plot colored shear flow overlay on airfoil in x/c, z/c (following airfoil)
%
% IMPORTANT: This is BASIC q' (boom method). For a fully closed multi-cell
% flexural shear flow, you'd add constant cell flows q0_i via compatibility.
% (We can add that after your plots are stable.)


%% ========== USER INPUTS ==========
csvFile = "ClarkY";      % <-- include extension
c = 160.1;                   % chord [mm]
Vz = V(1);                     % shear [N]

xF = 0.25;                   % front spar x/c
xR = 0.75;                   % rear spar x/c

% Skin/web thickness (only needed later if you want tau=q/t)
t_skin = 0.038;              % monokote [mm]
t_web  = 1.0;                % balsa/ply web [mm]

% Carbon tube assumption
Do = 8;                      % [mm]
Di = 6;                      % [mm]
A_tube = (pi/4)*(Do^2 - Di^2);

% Balsa stringer sizing rule
b_str = 3.0;                 % [mm]
k_str = 0.30;
hmax_str = 6.0;

% Stringer x/c positions
x_str_top = [0.125 0.35 0.50 0.65];
x_str_bot = [0.125 0.4167 0.5833];

% Plot controls
showLabels = true;
Nsamp = 120;

%% ========== 1) READ AIRFOIL ==========
M = readmatrix(csvFile);
x_all = M(:,1); z_all = M(:,2);
mask = ~(isnan(x_all) | isnan(z_all));
x_all = x_all(mask); z_all = z_all(mask);

% LE index
[~, iLE] = min(x_all);

% split
x1 = x_all(1:iLE);   z1 = z_all(1:iLE);
x2 = x_all(iLE:end); z2 = z_all(iLE:end);

% decide upper/lower by mean z
if mean(z1) > mean(z2)
    xU = x1; zU = z1;  xL = x2; zL = z2;
else
    xU = x2; zU = z2;  xL = x1; zL = z1;
end

% unique + sort
[xU, iu] = unique(xU, 'stable'); zU = zU(iu);
[xL, il] = unique(xL, 'stable'); zL = zL(il);
[xU, su] = sort(xU); zU = zU(su);
[xL, sl] = sort(xL); zL = zL(sl);

% interpolants (z/c)
zu_c = @(xc) interp1(xU, zU, xc, 'linear', 'extrap');
zl_c = @(xc) interp1(xL, zL, xc, 'linear', 'extrap');

%% ========== 2) BUILD BOOMS ==========
t_avail = @(xc) (zu_c(xc) - zl_c(xc))*c;                % [mm]
A_balsa = @(xc) b_str * min(k_str*t_avail(xc), hmax_str);

booms = struct('name',{},'xc',{},'x_mm',{},'z_mm',{},'B',{});
addBoom = @(name,xc,zmm,B) struct('name',string(name),'xc',xc,'x_mm',xc*c,'z_mm',zmm,'B',B);

% spars (top/bot)
booms(end+1) = addBoom("F_spar_top", xF, zu_c(xF)*c, A_tube);
booms(end+1) = addBoom("F_spar_bot", xF, zl_c(xF)*c, A_tube);
booms(end+1) = addBoom("R_spar_top", xR, zu_c(xR)*c, A_tube);
booms(end+1) = addBoom("R_spar_bot", xR, zl_c(xR)*c, A_tube);

% stringers
for k = 1:numel(x_str_top)
    xc = x_str_top(k);
    booms(end+1) = addBoom("StrTop_"+k, xc, zu_c(xc)*c, A_balsa(xc));
end
for k = 1:numel(x_str_bot)
    xc = x_str_bot(k);
    booms(end+1) = addBoom("StrBot_"+k, xc, zl_c(xc)*c, A_balsa(xc));
end

%% ========== 3) Iy (boom idealization) ==========
Bvec = [booms.B]';
zvec = [booms.z_mm]';
zbar = sum(Bvec .* zvec) / sum(Bvec);
zc   = zvec - zbar;
Iy   = sum(Bvec .* (zc.^2));

fprintf("Iy (boom idealization) = %.6e mm^4\n", Iy);

%% ========== 4) STATIONS PER CELL ==========
stations1 = buildStationsCell1(booms, xF, zu_c, zl_c, c, zbar);
stations2 = buildStationsCell2(booms, xF, xR, zu_c, zl_c, c, zbar);

%% ========== 5) BASIC q' (boom method) ==========
[qPanels1, ~] = shearFlowPanelsFromStations(stations1, Vz, Iy);
[qPanels2, ~] = shearFlowPanelsFromStations(stations2, Vz, Iy);

qAll = [abs(qPanels1(:)); abs(qPanels2(:))];
qMin = min(qAll);
qMax = max(qAll);

%% ========== 6) PLOT (x/c, z/c) ==========
figure(); 
hold on; grid on; axis equal;
title(sprintf("Flexural Shear Flow Overlay (Vz = %.1f N)", Vz));
xlabel("x/c"); ylabel("z/c");

% Outline
xc_plot = linspace(0,1,600);
plot(xc_plot, zu_c(xc_plot), 'k-', 'LineWidth', 1.2);
plot(xc_plot, zl_c(xc_plot), 'k-', 'LineWidth', 1.2);

% Webs
plot([xF xF], [zu_c(xF) zl_c(xF)], 'k-', 'LineWidth', 2.0);
plot([xR xR], [zu_c(xR) zl_c(xR)], 'k-', 'LineWidth', 2.0);

% Booms
scatter([booms.xc], [booms.z_mm]/c, 35, 'filled');
for i=1:numel(booms)
    text(booms(i).xc, booms(i).z_mm/c, "  "+booms(i).name, 'FontSize', 8);
end

% --- IMPORTANT: set colormap + force color limits BEFORE drawing ---
colormap(parula(256));
set(gca,'CLim',[qMin qMax]);

% Dummy mappable object so colorbar ALWAYS appears
hDummy = scatter(nan,nan,1,nan,'filled');
set(hDummy,'Visible','off');

% Draw both cells (colored, follows airfoil)
plotColoredPanelsAirfoil_xc(stations1, qPanels1, qMin, qMax, showLabels, c, zu_c, zl_c, Nsamp);
plotColoredPanelsAirfoil_xc(stations2, qPanels2, qMin, qMax, showLabels, c, zu_c, zl_c, Nsamp);

cb = colorbar;
cb.Label.String = "|q'| (N/mm)";

fprintf("\nDone. Plot shows BASIC flexural shear flow q'(s).\n");

%% ===================== LOCAL FUNCTIONS =====================

function stations = buildStationsCell1(booms, xF, zu_c, zl_c, c, zbar)
% LE cell:
% front spar TOP -> upper to LE -> lower back -> front spar BOT -> web up

top = booms([booms.z_mm] > 0 & [booms.xc] <= xF+1e-9);
bot = booms([booms.z_mm] < 0 & [booms.xc] <= xF+1e-9);

% Remove any accidental booms at x/c~0 (prevents weird LE behavior)
tol = 1e-6;
top = top([top.xc] > tol);
bot = bot([bot.xc] > tol);

[~,idx] = sort([top.xc],'descend'); top = top(idx);  % xF -> 0
[~,idx] = sort([bot.xc],'ascend');  bot = bot(idx);  % 0 -> xF

stations = struct('name',{},'x_mm',{},'z_mm',{},'B',{},'zc',{},'isSkin',{},'side',{});

for i=1:numel(top)
    stations(end+1) = makeStationFromBoom(top(i), "top", zbar);
end

% Two LE stations at SAME point so no fake segment appears
% (this is the key "nose fix")
stations(end+1) = struct('name',"LE_top", 'x_mm',0, 'z_mm',0, 'B',0, 'zc',0, ...
                         'isSkin', true, 'side',"top");
stations(end+1) = struct('name',"LE_bot", 'x_mm',0, 'z_mm',0, 'B',0, 'zc',0, ...
                         'isSkin', true, 'side',"bot");

for i=1:numel(bot)
    stations(end+1) = makeStationFromBoom(bot(i), "bot", zbar);
end

% Close via front web up to top
stations(end+1) = struct('name',"F_web_top", 'x_mm',xF*c, 'z_mm',zu_c(xF)*c, 'B',0, 'zc',0, ...
                         'isSkin', false, 'side',"web");
end

function stations = buildStationsCell2(booms, xF, xR, zu_c, zl_c, c, zbar)
% Middle cell:
% front spar TOP -> upper to rear spar TOP -> rear web down -> lower back -> front web up

top = booms([booms.z_mm] > 0 & [booms.xc] >= xF-1e-9 & [booms.xc] <= xR+1e-9);
bot = booms([booms.z_mm] < 0 & [booms.xc] >= xF-1e-9 & [booms.xc] <= xR+1e-9);

[~,idx] = sort([top.xc],'ascend');  top = top(idx);     % xF -> xR
[~,idx] = sort([bot.xc],'descend'); bot = bot(idx);     % xR -> xF

stations = struct('name',{},'x_mm',{},'z_mm',{},'B',{},'zc',{},'isSkin',{},'side',{});

for i=1:numel(top)
    stations(end+1) = makeStationFromBoom(top(i), "top", zbar);
end

% Rear web down
stations(end+1) = struct('name',"R_web_bot", 'x_mm',xR*c, 'z_mm',zl_c(xR)*c, 'B',0, 'zc',0, ...
                         'isSkin', false, 'side',"web");

for i=1:numel(bot)
    stations(end+1) = makeStationFromBoom(bot(i), "bot", zbar);
end

% Front web up
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

q = 0;
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

    % Skip zero-length segments
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

    % Draw as a "surface edge" so MATLAB does NOT close the polygon (no double lines)
    X = [xcs; xcs];
    Y = [zcs; zcs];
    Z = zeros(size(X));
    C = qv * ones(size(X));

    h = surface(X, Y, Z, C, ...
        'FaceColor','none', ...
        'EdgeColor','interp', ...
        'LineWidth',4);

    set(h,'CDataMapping','scaled');

    if showLabels
        text(mean(xcs), mean(zcs), sprintf(" %.2f", qPanels(k)), 'FontSize', 8);
    end
end

% Keep global color scaling consistent
set(gca,'CLim',[qMin qMax]);
end

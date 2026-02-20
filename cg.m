% mass in kg, x cg in m
balsa_density = 130; % kg/m^3
carbon_density = 2 *1000; % g / cm^3


%% Fuselage

% 124.1 mm = front of fuselage, 521 mm back of fuselage
x_fuselage_start = .1241;
x_fuselage_back = .521;
n_hoops = 10;
a_hoops = 5233 / 1000000; % 5233 mm^2
t_hoops = .002; % 2 mm 

a_fuselage_l1 = 9289 / 1000000; % there are 4 of these
a_fuselage_l2 = 10828 / 1000000; % there are 2 of these
a_fuselage_l3 = 11864 / 1000000; % there are 2 of these

m_fuselage_hoops = n_hoops * a_hoops * t_hoops * balsa_density;
m_fuselage_lengthwise = t_hoops * balsa_density * (a_fuselage_l1 + ...
    a_fuselage_l2 + a_fuselage_l3);
m_fuselage = m_fuselage_hoops + m_fuselage_lengthwise;
x_cg_fuselage = .5 * (x_fuselage_start + x_fuselage_back); % 300 mm

%% Prop System
x_cg_prop = 0.022;
m_prop = .023;
x_cg_motor = .05; % 50 mm
m_motor = .119;
x_cg_ESC = x_fuselage_start + .005;
l_ESC = 0.073;
m_ESC = .065;

x_cg_battery = x_cg_ESC + l_ESC;
l_battery = .14;
m_battery = .3777;

x_cg_Pixhawk = x_cg_ESC + l_battery; 
m_Pixhawk = .035;

m_prop_system = m_prop + m_motor + m_ESC + m_Pixhawk + m_battery;
x_cg_prop_system = (x_cg_prop*m_prop + x_cg_motor * m_motor + ...
    x_cg_ESC * m_ESC + x_cg_Pixhawk * m_Pixhawk + ...
    x_cg_battery * m_battery) / m_prop_system;

%% Wing
wing_spar_area = 831 / 1000000; % m^2
wing_spar_thickness = .002; % mm
wing_spar_volume = wing_spar_thickness * wing_spar_area; % mm^2 area

m_ribs_w = wing_spar_volume * balsa_density;

m_spar_w = (37.6991 / 1000000) * b * carbon_density;
n_ribs_w = 15;

x_cg_top = .2608; % spinner to quarter chord
m_top = m_spar_w + m_ribs_w * n_ribs_w;
x_cg_bottom = .32412; % spinner to quarter chord
m_bottom = m_top;

m_wings = m_top + m_bottom;
x_cg_wings = (x_cg_top * m_top + x_cg_bottom * m_bottom) / m_wings;

% %% Tail
% n_ribs_vt = ;
% n_ribs_ht = ;
% 
% m_rib_vt = ;
% m_rib_ht = ;
% 
% x_cg_vt = ;
% m_vt = n_ribs_vt * m_rib_vt;
% x_cg_ht = ;
% m_ht = n_ribs_ht + m_rib_ht;
% x_cg_tail_boom = number + lt / 2;
% m_tail_boom = 
% 
% m_tail = m_vt + m_ht + m_tail_boom;
% x_cg_tail = (x_cg_vt * m_vt + x_cg_ht * m_ht + x_cg_tail_boom * m_tail_boom) / m_tail;

m_tail = Sh_S * m_wings;
x_cg_tail = 1.3012;

% %% Landing Gear
x_cg_noselg = .06 + x_fuselage_start;
m_nose_lg = .166; % kg
x_cg_mainlg = .35 + x_fuselage_start;
m_main_lg = .391; % kg

m_lg = m_main_lg + m_nose_lg;
x_cg_lg = (x_cg_noselg * m_nose_lg + x_cg_mainlg * m_main_lg) / m_lg;

m_plane = m_lg + m_tail + m_wings + m_prop_system + m_fuselage;
x_cg = (m_lg * x_cg_lg + m_tail * x_cg_tail + m_wings * x_cg_wings + ...
    m_prop_system * x_cg_prop_system+ m_fuselage * x_cg_fuselage) / m_plane;
% calculation for the location of the landing gear

alpha_TO = (CL_R - CL_0)/(CL_a);
alpha_TO_deg = rad2deg(alpha_TO);

alpha_TB_deg = alpha_TO_deg + 5; % add 5 degrees to the TO alpha for tipback angle
alpha_TB = deg2rad(alpha_TB_deg); % convert tipback angle to radians

x_mg = (x_cg_chosen) + (0.21*tan(alpha_TB)); % the x_cg needs to be fixed!

nose_load = 0.2; % 20% of the wieght should be on the nose gear

x_nose = x_mg - (x_mg - x_cg_chosen)/nose_load;
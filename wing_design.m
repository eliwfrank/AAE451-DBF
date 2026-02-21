fprintf("\nWING SIZING ------------------------\n")

wing_area_total = weight_TO / W_S;

wing_area_ind = wing_area_total / 2; % because of the biplane we need individual

b = sqrt(AR_wing * wing_area_ind); % AR = b^2/A where b is span
c = wing_area_ind / b; % chord

wing_thickness_cm = 100* 0.117*c; % thickness is 11.7% of chord for the clark y

fprintf("Wing Chord length: %.4f [m]\nWing Span: %.4f [m]\n",c,b)
fprintf("Wing thickness: %.4f [m]\n",wing_thickness_cm);
fprintf("Wing Area: %.4f [m^2]\n",wing_area_total);
wing_area_total = W_int / W_S;

wing_area_ind = wing_area_total / 2; % because of the biplane we need individual

b = sqrt(AR * wing_area_ind); % AR = b^2/A where b is span
c = wing_area_ind / b; % chord

fprintf("Chord length: %.4f [m]\nSpan: %.4f [m]\n",c,b)
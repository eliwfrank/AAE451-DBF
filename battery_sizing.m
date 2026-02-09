%% Battery Sizing

S = 4; %Battery Series               
I_avg = 30; %Average Current Draw - this is important to get right, could break up into current draw at different flight stages      
flight_time_hr = (t_TO + t_M + t_CR + t_CL) / (60*60); %flight time

usable_fraction = 0.80;   % Only use 80% of capacity (LiPo safety)
safety_margin = 1.5;     % 10% extra capacity margin

Ah_required = I_avg * flight_time_hr;   
Ah_required = Ah_required / usable_fraction; 
Ah_required = Ah_required * safety_margin; 

mAh_required = Ah_required * 1000;

fprintf("Battery selection for %dS system:\n",S)
fprintf("Required capacity ≈ %.0f mAh\n",mAh_required)

V_nom = 3.6 * S;  % nominal LiPo voltage - ask Tom if we can go up to 4.2 charge voltage?
P_avg = V_nom * I_avg;

fprintf("Estimated cruise power ≈ %.1f W\n",P_avg)


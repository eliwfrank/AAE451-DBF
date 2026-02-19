%% Battery Sizing
fprintf("\nBATTERY SIZING ------------------------\n")
S = 5; %Battery Series 
voltage = 3.7;
Pavg = 550; %average flight power in watts
I_avg = Pavg/(S*voltage); %Average Current Draw - this is important to get right, could break up into current draw at different flight stages      
flight_time_hr = (t_TO + t_M + t_CR + t_CL) / (60*60); %flight time

usable_fraction = 0.80;   % Only use 80% of capacity (LiPo safety)
safety_margin = 1.5;     %1.5x flight time

Ah_required = I_avg * flight_time_hr;   
Ah_required = Ah_required / usable_fraction; 
Ah_required = Ah_required * safety_margin; 

mAh_required = Ah_required * 1000;

fprintf("Battery selection for %dS system:\n",S)
fprintf("Required capacity ≈ %.0f mAh\n",mAh_required)





% Intertidal seagrass model that simulate the changes of seagrass growth rate when response to periodic tidal inundation and exposure
% Main code for baseline scenario

clc
clear all
%% Simulate the intertidal cycle of air exposure and submersion

% Calculate the water depth, m
A_M2 = 0.6;  % amplitude of M2 tidal constituent，m
T_M2 = 0.5175; % period of M2 tidal constituent, d
A_S2 = 0.2;  % amplitude of S2 tidal constituent，m
T_S2 = 0.5; %  period of S2 tidal constituent, d
Z_b_values = [-3.0, -2.9, -2.8, -2.7, -2.6, -2.5, -2.4,-2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1,-1.0,-0.9,-0.8, -0.7, -0.6,-0.5,-0.4,-0.3,-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0]; % meadow elevation relative to mean sea level, m

for j = 1:length(Z_b_values) % Loop through each Z_b value
    Z_b = Z_b_values(j);
    
T = 15; % model simulation time, d
delta_t = 1/1440; % model simulation interval, d
t = linspace(0, T, T/delta_t + 1); % time array in days

Z_w = A_M2*cos(2*pi*t/T_M2) + A_S2*cos(2*pi*t/T_S2); % calculate the water level, m

for i = 1:length(t) % calculate the water depth, m
    if Z_w(i)>Z_b
    Z_t(i)= Z_w(i)-Z_b;
    else
    Z_t(i)=0;  
    end
end

% Calculate the air temperature, ℃
T_ave = 25;  % mean daily air temperature,℃
T_var = 5; % daily air temperature from its mean value, ℃
t_solar = 0.083; % time between solar noon and when the maximum temperature occurs in the afternoon in days, d 

T_t = T_ave-T_var*cos(2*pi*(t-t_solar)); % air temperature at time t, ℃


% Calculate the benthic light irradiation
for i = 1:T+1 % calculate daily t_rise and t_set
    t_rise(i) = 0.2708 + (i-1); % sunrise time in days,d
    L=0.5; % total daylight in days, d
    t_set(i) = t_rise(i) + L; % sunset time in days,d
end

I_0 = 60; % daily average water surface light, mol m-2 d-1
K_d = 0.5; % light attenuation coefficient in the water column

for i = 1:length(t) % dalily average benthic light irradiance, mol m-2 d-1
    if Z_t(i) > 0
      I_ave(i) = I_0*exp(-K_d*Z_t(i));
    else
      I_ave(i) = I_0;
    end
end

for i = 1:length(t) % benthic light irradiance at time t, mol m-2 d-1
    t_rise_rep = reshape(repmat(t_rise, 1/delta_t, 1), 1, []); % repeat sunrise time every 24 hours
    t_set_rep = reshape(repmat(t_set, 1/delta_t, 1), 1, []); % repeat sunrise time every 24 hours

    if t(i) >= t_rise_rep(i) && t(i) <= t_set_rep(i) % time during daylight 
        I_t(i) =(pi/(2*L))*I_ave(i)*sin(pi/L*(t(i)-t_rise_rep(i))); % calculate the light irradiance using a sine function
    else
        I_t(i) = 0;
    end
end

% Calculate the Desiccation coefficient for seagrass
for i = 1:length(t)
    k(i) = 0.451*(T_t(i)-20)+11.617; % desiccation coefficient at temperature T, d-1
end

% Calculate air exposure duration
t_air = zeros(size(t)); 
for i = 2:length(t)
    if Z_t(i-1) > 0 && Z_t(i) == 0
     t_0 = interp1(Z_t(i-1:i), t(i-1:i), 0); % calculate the t_0 that makes water depth Z_t=0
     t_air(i) = t(i) - t_0;   
    elseif Z_t(i-1) > 0 && Z_t(i) > 0
        t_air(i) = 0;
    elseif Z_t(i-1) == 0 && Z_t(i) == 0
     t_air(i) = t_air(i-1) + t(i)-t(i-1);  % calculate air exposure duration at each time interval, d 
    end
end

% Calculate Relative water content of seagrass leaves at each time interval
for i = 1:length(t)
    if Z_t(i) > 0
        RWC_t(i) = 1;
    else
        RWC_t(i) = 1*exp(-k(i)* t_air(i));
    end
end

%% Calculate intertidal seagrass growth rate and environmental stress (desiccation stress, light deprivation stress and overall stress)

% Define parameters for seagrass growth function (sigmoidal curve model of f_RWC (RWC))
u_max = 0.04; % maximum seagrass growth rate, d-1
a =0.01; % efficiency of light utilisation for seagrass growth at low light, d-1/(mol m-2 d-1)
R = 15.043; % shape parameter in sigmoidal curve model of f_RWC (RWC) 
RWC_h = 0.381; % the value of RWC that attains half of the effective quantum yield in sigmoidal curve model of f_RWC (RWC)

for i = 1:length(t)
   fRWC(i) = 1/(1+exp(-R*(RWC_t(i)-RWC_h))); % calculate the f_RWC(RWC) 
   fI(i) = a * I_t(i)/ (a * I_t(i) + u_max); % calculate the f_I(I) 
   u1(i) = u_max * fI(i) * fRWC(i);  % calculate the seagrass growth rate dependent on I,RWC, multiplicative formulation
   u2(i) = u_max * min(fI(i), fRWC(i)); % calculate the seagrass growth rate dependent on I,RWC, law of minimum formulation
   u3(i) = u_max * fI(i); % calculate the growth rate dependent on I

 mean_fRWC(j) = 1-mean(fRWC, "all"); % average desiccation stress at each meadow elevation during simulation
 mean_fI(j) = 1-mean(fI, "all"); % average light deprivation stress at each meadow elevation during simulation 
 mean_overall1(j) = 1- mean(u1, "all")/u_max % overall stress with multiplicative formulation at each meadow elevation during simulation
 mean_overall2(j) = 1- mean(u2, "all")/u_max % overall stress with law of minimum formulation at each meadow elevation during simulation


 % Store the meadow elevation Z_b, desiccation stress mean_fRWC, light deprivation stress mean_fI and overall stress values mean_overall1/2 values in the output array
    output_data(1, j) = Z_b;
    output_data(2, j) = mean_fRWC(j);
    output_data(3, j) = mean_fI(j);
    output_data(4, j) = mean_overall1(j);
    output_data(5, j) = mean_overall2(j);
end
end

%% BALLOOON LAB
clear;
clc;
close all;

%% Variables

% General
T_atm = 295.15;               %K
P_atm = 84800;                %Pa
R_air = 287;                  %J/kgK
g = 9.81;                     %m/s^2

% Helium Balloon
R_He = 207.69;                %J/kg*K
m_helium_bal = 0.009569;      %kg
m_helium_payload = 0.001086;  %kg
displace_tank = .098;         %m
length_tank = 0.495;          %m
width_tank = 0.248;           %m


% Hot Air Balloon
V_hot_air = 0.1313902; %m^3
m_hot_air_shell = .033; %kg
T_max = 522; %K

%% Helium Calculations
% Calculate the mass of the Helium Ideal Gas Law
V_helium_bal = displace_tank * length_tank * width_tank;   %m^3
density_air = P_atm / (R_air * T_atm);                     %m^3/kg
mass_He = density_air*V_helium_bal - m_helium_bal - m_helium_payload;  %kg

%disp(V_helium_bal)
fprintf("Helium Volume (m^3): %.4f \n", V_helium_bal);
fprintf("Air Density (m^3/kg): %.4f \n", density_air);
fprintf("Helium Mass (kg): %.4f \n", mass_He);

mass_total = mass_He + m_helium_payload + m_helium_bal; %kg

%Finding the mass of helium through the Force Balance Equation
% 0=ForceBuoyant - WeightBalloon - WeightPayload - WeigthHelium
ForceBouyant = 9.81*(m_helium_bal + m_helium_payload + mass_He); %N

%Plotting the Volume over Altitude (0-30000) = [1:3000];
%Creating a Vector for the height traveled
data_vol = 1:30000;
T_vec = 1:30000;
for i = data_vol
    [T,a,P,rho] = atmoscoesa(i);
    T_vec(i) = T;
    data_vol(:,i) = balloon_mass_func(mass_total,T,P,R_He);
end

%% Hot Air Calculations

T_hotair = findInternalTemp(P_atm, T_atm, V_hot_air, m_hot_air_shell, R_air);

fprintf("Hot Air Balloon Temperature (K): %.4f \n", T_hotair);

T_hotair_vec = 1:10000;
for i = T_hotair_vec
    [Temp, vSound, Pressure, Density] = atmoscoesa(i);
    T_hotair_vec(:,i) = findInternalTemp(Pressure, Temp, V_hot_air, m_hot_air_shell, R_air);
end

T_hotair_vec(T_hotair_vec(:) > T_max) = NaN;

%% Graphing Hot Air Data

figure();
plot(1:10000,T_hotair_vec, "LineWidth", 1); hold on;
yline(max(T_hotair_vec), "--", "Max Temperature", "LabelHorizontalAlignment","left");
xline(find(T_hotair_vec == max(T_hotair_vec)), "--", "Max Altitude","LabelVerticalAlignment","bottom");
title("Temperature of the Hot Air Balloon vs Altitude")
xlabel("Altitude (m)")
ylabel("Temperature (\circ K)")
hold off;

%% Graphing Volume Data

%Graphing the Volume vs Altitude
Temp_max = max(T_vec);
Volume = figure();
plot(1:30000,data_vol,'blue', 'LineWidth',1); hold on;
line([0 find(T_hotair_vec == max(T_hotair_vec))], [V_hot_air V_hot_air]);
xline(find(T_hotair_vec == max(T_hotair_vec)), "--");
title('System Volumes vs Altitude')
ylabel('Volume (m^3)')
xlabel('Altitude (m)')
legend("Helium Balloon Volume","Hot Air Balloon Volume", "Max Hot Air Balloon Altitude")
hold off;

%% Funtions

% This balloon mass function calculates the mass using the ideal gas law
% rho = mass*pressure/(R*T)
function [Volume] = balloon_mass_func(Mass,Temperature, Pressure, R_Constant)
    Volume = R_Constant*Temperature*Mass/Pressure;
end

% Calculates internal temp of the hot air balloon using ideal gas law and
% the neutral buoyancy equation
function [internalTemp] = findInternalTemp(P_atm, T_atm, V_hot_air, m_hot_air_shell, R_air)
    internalTemp = (P_atm * T_atm * V_hot_air)/(-m_hot_air_shell * R_air * T_atm + P_atm * V_hot_air);
end
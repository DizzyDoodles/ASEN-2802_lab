%% BALLOOON LAB
clear;
clc;
close all;

%% Helium Balloon Data
m_bal = 0.009569;      %Kilograms
m_payload = 0.001086;  %Kilograms
Patm = 84000;          %Pa
Tatm = 295;            %degrees K
Displace_Tank = .098;  %meters
LengthTank = 0.495;    %meters
WidthTank = 0.248;     %meters
g = 9.81;              %m/s^2
R_He = 207.69;         %J/kg*K
R_air = 287;           %J/kg*K

%% Calculations
% Calculate the mass of the Helium Ideal Gas Law
VolBal = Displace_Tank*LengthTank*WidthTank; %m^3
density_air = Patm/(R_air * Tatm); %m^3/kg
mass_He = density_air*VolBal - m_bal - m_payload; %kg

disp(VolBal)
disp(density_air)
disp(mass_He)

mass_total = mass_He + m_payload + m_bal; %Kilograms


%Finding the mass of helium through the Force Balance Equation
% 0=ForceBuoyant - WeightBalloon - WeightPayload - WeigthHelium
ForceBouyant = 9.81*(m_bal + m_payload+mass_He); %Newtons

%Plotting the Volume over Altitude (0-30000) = [1:3000];
%Creating a Vector for the height traveled
data_vol = 1:30000;
T_vec = 1:30000;
for i = data_vol
    [T,a,P,rho] = atmoscoesa(i);
    T_vec(i) = T;
    data_vol(:,i) = balloon_mass_func(mass_total,T,P,R_He);
end

%% Graphing the Data
%Graphing the Volume vs Altitude
Temp_max = max(T_vec);
Volume = figure();
plot(1:30000,data_vol,'blue', 'LineWidth',1)
title('Volume vs Altitude')
ylabel('Volume(m^3')
xlabel('Altitude (m)')
legend('The Legendary Line', 'location', 'southeast')
hold off

%Graphing the Temperature vs Altitude
Temperature = figure();
plot(T_vec,1:30000, 'red', 'LineWidth', 1)
title('Temperature vs Altitude')
ylabel('Altitude (m)')
xlabel('Temperature (k)')
legend('The Legendary Line 2')



%% Functions
%This balloon mass function calculates the mass using the ideal gas law
%rho = mass*pressure/(R*T)
function [Volume] = balloon_mass_func(Mass,Temperature, Pressure, R_Constant)
    Volume = R_Constant*Temperature*Mass/Pressure;
end



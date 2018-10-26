%% AE 4451 Jet & Rocket Propulsion Project
%10/26/2018
%Authors: Loren Isakson, Brandon Campanile, Matthew Yates
clear;clc
close all
%% Inputs
Ta = 0;    %ambient temperature, Kelvin
Pa = 0;    %ambient pressure, Pa
Pf = 0;    %fuel storage pressure
M = 0;      %flight mach number
Prc = 0;   %compressor stagnation pressure ratio
Prf = 0;   %fan stagnation pressure ratio
f = 0;      %main burner fuel-air ratio
fab = 0;   %afterburner fuel-air ratio
%% Inputs
T_a = 0;    %ambient temperature
P_a = 0;    %ambient pressure
P_f = 0;    %fuel storage pressure
M = 0;      %flight mach number
Pr_c = 0;   %compressor stagnation pressure ratio
Pr_f = 0;   %fan stagnation pressure ratio
f = 0;      %main burner fuel-air ratio
f_ab = 0;   %afterburner fuel-air ratio
beta = 0;   %bypass ratio
b = 0;      %bleed ratio

global R   %global variables
R = 8314/28.8;  %universal gas constant / molecular weight of the gas

%% Function Library

%Ambient conditions provided T_a and P_a

function [To1, Po1] = diffuser(Ta, Pa, M)
%static ambient temp, press, mach number, gamma, adiabatic efficiency
gamma = 1.4;
nd = 0.92;
To1 = Ta.*(1+0.5.*(gamma-1).*M.^2);
Po1 = Pa.*(1+nd.*(To1/Ta - 1)).^(gamma/(gamma-1));
end

function [To2, Po2, wf_ma] = fan(To1, Po1, Prf, beta)
%stagnation temp, press, pressure ratio, polytropic efficiency
gamma = 1.4;
npf = 0.9;
Cp = gamma*R/(gamma-1);
if Prf < 1.1 || Prf > 1.5
    error('Fan pressure ratio is constrained 1.1 <= Pr <= 1.5');
end 
To2 = To1.*(Pr).^((gamma-1)/gamma/npf);
Po2 = Po1.*Pr;
wf_ma = Cp*(1+beta)*(To2-To1); %work done by the fan on the fluid
end

function [To3, Po3, wc_ma] = compressor(To2, Po2, Prc, Prf)
gamma = 1.38;
npc = 0.9;
Cp = gamma*R/(gamma-1);
if Prc < 60/Prf
    error('Compressor pressure ratio must be less than 60/fan pressure ratio');
end
Po3 = Po2*Prc;
To3 = To2*(Prc).^((gamma-1)/gamma/npc);
wc_ma = Cp*(To3-To2);
end

function [To4, Po4] = burner(To3, Po3, Tmax, deltaH)
gamma = 1.33;
nb = 0.99;
Prb = 0.98;
Cp = gamma*R/(gamma-1);

end

function [To5_1, Po5_1] = turbine(To4, Po4)

end

function [To5_m, Po5_m] = turbineMixer(To5, Po5)

end

function [To5_2, Po5_2] = fanTurbine(To5_m, Po5_m)

end

function [To6, Po6] = afterburner(To5_2, Po5_2)

end

function [Toe, Poe] = coreNozzle(To6, Po6)

end

function [Toef, Poef] = fanNozzle(Toe, Poe)

end

function [To7, Po7] = nozzleMixer(Toef, Poef)

end

function [Toec, Poec] = combinedNozzle(To7, Po7)

end

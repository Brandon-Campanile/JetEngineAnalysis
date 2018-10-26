%% AE 4451 Jet & Rocket Propulsion Project
%10/26/2018
%Authors: Loren Isakson, Brandon Campanile, Matthew Yates
clear;clc
%% Inputs
T_a = 0;    %ambient temperature
P_a = 0;    %ambient pressure
P_f = 0;    %fuel storage pressure
M = 0;      %flight mach number
Pr_c = 0;   %compressor stagnation pressure ratio
Pr_f = 0;   %fan stagnation pressure ratio
f = 0;      %main burner fuel-air ratio
f_ab = 0;   %afterburner fuel-air ratio
Beta = 0;   %bypass ratio
b = 0;      %bleed ratio

%% Function Library

%Ambient conditions provided T_a and P_a
function table = engineAnalysis(Ta, Pa, Pf, M, Prc, Prf, B, b, f, fab)
% Ta = ambient temperature [Kelvin]; Pa = ambient pressure [kPa]; Pf = fuel storage pressure [kPa]; M = flight mach number
% Prc = compressor stagnation pressure ratio; Prf = fan stagnation; pressure; ratio; B = bypass ratio; b = bleed ratio
% f = main burner fuel-air ratio; fab = afterburner fuel-air ratio
if b = 0;
    


end


function [To1, Po1] = diffuser(Ta, Pa, M)
%static ambient temp, press, mach number, gamma, adiabatic efficiency
gamma = 1.4;
nd = 0.92;
To1 = Ta.*(1+0.5.*(gamma-1).*M.^2);
Po1 = Pa.*(1+nd.*(To1/Ta - 1)).^(gamma/(gamma-1));
end

function [To2, Po2] = fan(To1, Po1, Pr)
%stagnation temp, press, pressure ratio, polytropic efficiency
gamma = 1.4;
np = 0.9;
MW = 28.8;
R = 8314/MW;
Cp = gamma*R/(gamma-1);
if Pr < 1.1 || Pr > 1.5
    error('Fan pressure ratio is constrained 1.1 <= Pr <= 1.5');
end 
To2 = To1.*(Pr).^((gamma-1)/gamma/np);
Po2 = Po1.*Pr;
end

function [To3, Po3] = compressor(To2, Po2, gamma, n, MW)

end

function [To4, Po4] = burner(To3, Po3, deltaH)

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
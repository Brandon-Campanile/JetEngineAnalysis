%% AE 4451 Jet & Rocket Propulsion Project
%10/26/2018
%Authors: Loren Isakson, Brandon Campanile, Matthew Yates

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

function [T_o1, P_o1] = inlet(M, T_a, P_a, gamma)
    T_o1 = T_a.*(1+0.5.*(gamma-1).*M.^2);
    P_o1 = P_a.*(1+0.5.*(gamma-1).*M.^2).^(gamma/(gamma-1));
end

function [T_o2, P_o2] = diffuser(T_o1, P_o1, gamma, n, MW)

end

function [T_o3, P_o3] = compressor(T_o2, P_o2, gamma, n, MW)

end

function [T_o4, P_o4] = combustor(T_o3, P_o3







%
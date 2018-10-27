%% AE 4451 Jet & Rocket Propulsion Project
%10/26/2018
%Authors: Loren Isakson, Brandon Campanile, Matthew Yates
clear;clc
close all
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

global R_   %global variables
R_ = 8314;  %universal gas constant

%% Component Functions

%Ambient conditions provided T_a and P_a
function table = engineAnalysis(Ta, Pa, Pf, M, Prc, Prf, B, b, f, fab, ab, mix)
% Ta = ambient temperature [Kelvin]; Pa = ambient pressure [kPa]; Pf = fuel storage pressure [kPa]; M = flight mach number
% Prc = compressor stagnation pressure ratio; Prf = fan stagnation; pressure; ratio; B = bypass ratio; b = bleed ratio
% f = main burner fuel-air ratio; fab = afterburner fuel-air ratio; ab = afterburner boolean; mix = nozzle mixing boolean
if Prf
    [To1, Po1] = diffuser(Ta, Pa, M);
    [To2, Po2, wf_ma] = fan(To1, Po1, Prf, B);
    [To3, Po3, wc_ma] = compressor(To2, Po2, Prc, Prf);
    [To4, Po4, wp_ma] = burner(Po3, b, f);
    [To5_1, Po5_1] = turbine(To4, Po4, wc_ma, wp_ma, b, f);
    [To5_m, Po5_m] = turbineMixer(To5, Po5);
    [To5_2, Po5_2] = fanTurbine(To5_m, Po5_m);
    if ab && mix
        [To6, Po6] = afterburner(Po5_2, Tmax_ab, Prab);
        [To7, Po7] = nozzleMixer(Toef, Poef);
        [Toec, Poec] = combinedNozzle(To7, Po7);
    elseif ab && ~mix
        [To6, Po6] = afterburner(Po5_2, Tmax_ab, Prab);
        [Toef, Poef] = fanNozzle(Toe, Poe);
        [Toe, Poe] = coreNozzle(To6, Po6);
    elseif ~ab && mix
        [To7, Po7] = nozzleMixer(Toef, Poef);
        [Toec, Poec] = combinedNozzle(To7, Po7);
    else
        [Toe, Poe] = coreNozzle(To6, Po6);
    end
else
    [To1, Po1] = diffuser(Ta, Pa, M);
    [To3, Po3, wc_ma] = compressor(To1, Po1, Prc, Prf);
    [To4, Po4] = burner(To3, Po3, Tmax, deltaH);
    [To5_1, Po5_1] = turbine(To4, Po4);
    [To5_m, Po5_m] = turbineMixer(To5_1, Po5_1);
    if ab
        [To6, Po6] = afterburner(Po5_2, Tmax_ab, Prab);
    else
        To6 = To5_m;
        Po6 = Po5_m;
    end
    [Toe, Poe] = coreNozzle(To6, Po6);
end
end

function [To1, Po1] = diffuser(Ta, Pa, M)
%static ambient temp, press, mach number, gamma, adiabatic efficiency
gamma = 1.4;
nd = 0.92;
To1 = Ta.*(1+0.5.*(gamma-1).*M.^2);
Po1 = Pa.*(1+nd.*(To1/Ta - 1)).^(gamma/(gamma-1));
end

function [To2, Po2, wf_ma] = fan(To1, Po1, Prf, beta, MW)
%stagnation temp, press, pressure ratio, polytropic efficiency
gamma = 1.4;
npf = 0.9;
Cp = gamma*(R_/MW)/(gamma-1);
if Prf < 1.1 || Prf > 1.5
    error('Fan pressure ratio is constrained 1.1 <= Pr <= 1.5');
end 
To2 = To1.*(Pr).^((gamma-1)/gamma/npf);
Po2 = Po1.*Pr;
wf_ma = Cp*(1+beta)*(To2-To1); %work done by the fan on the fluid
end

function [To3, Po3, wc_ma] = compressor(To2, Po2, Prc, Prf, MW)
gamma = 1.38;
npc = 0.9;
Cp = gamma*(R_/MW)/(gamma-1);
if Prc < 60/Prf
    error('Compressor pressure ratio must be less than 60/fan pressure ratio');
end
Po3 = Po2*Prc;
To3 = To2*(Prc).^((gamma-1)/gamma/npc);
wc_ma = Cp*(To3-To2);
end

function [To4, Po4, wp_ma] = burner(Po3, b, f, MW)
gamma = 1.33;
nb = 0.99;
Prb = 0.98;
Cp = gamma*(R_/MW)/(gamma-1);
Tomax = 1300; %kelvin
Cb = 700; %kelvin
bmax = 0.12;
To4 = Tomax + Cb*(b/bmax)^0.5;
Po4 = Po3*Prb;
%Fuel Pump
deltaP_inject = 550*10^3;
np = 0.35;
f_density = 780;
wp_ma = f*deltaP_inject/np/f_density;
end

function [To5_1, Po5_1] = turbine(To4, Po4, wc_ma, wp_ma, b, f, MW)
npt = 0.92;
gamma = 1.33;
Cp = gamma*(R_/MW)/(gamma-1);
To5_1 = To4 - (1/(Cp*(1+f-b))*(wc_ma + wp_ma));
TR = To5_1/To4;
Po5_1 = Po4*(1+((TR-1)^2)/(TR^(1/npt) - 1))^(gamma/(gamma-1));
end

function [To5_m, Po5_m] = turbineMixer(To5_1, Po5_1, To3, f, b, MW)
gamma = 1.34;
Cp = gamma*(R_/MW)/(gamma-1);
Po5_m = Po5_1;
To5_m = (b*To3 + (1 + f - b)*To5_1)/(1+f);
end

function [To5_2, Po5_2] = fanTurbine(To5_m, Po5_m, wf_ma, f, MW)
gamma = 1.33;
Cp = gamma*(R_/MW)/(gamma-1);
To5_2 = To5_m - wf_ma/(Cp*(1+f));
TR = To5_2/To5_m;
Po5_2 = Po5_m*(1 + ((TR-1)^2)/(TR-1))^(gamma/(gamma-1));
end

function [To6, Po6] = afterburner(Po5_2, Tmax_ab, Prab)
To6 = Tmax_ab;
Po6 = Prab*Po5_2;
end
<<<<<<< HEAD

function [Te, ue] = coreNozzle(To6, Po6, Pa)
=======
%%
function [Te, ue] = coreNozzle(To6, Po6, Pa, MW)
>>>>>>> c2925fbee04d1c1f09affd74e3c9d272a73fc902
gamma = 1.35;
nc = 0.95;
Cp = gamma*(R_/MW)/(gamma-1);
Te = To6*(1-nc*(1-(Pa/Po6)^((gamma-1)/gamma)));
ue = sqrt(2*Cp*(To6-Te));
end
<<<<<<< HEAD

function [Tef, uef] = fanNozzle(To2, Po2, Pa)
=======
%+
function [Tef, uef] = fanNozzle(To2, Po2, Pa, MW)
>>>>>>> c2925fbee04d1c1f09affd74e3c9d272a73fc902
gamma = 1.4;
nf = 0.97;
Cp = gamma*(R_/MW)/(gamma-1);
Tef = To2*(1-nf*(1-(Pa/Po2)^((gamma-1)/gamma)));
uef = sqrt(2*Cp*(To2-Tef));
end

function [To7, Po7] = nozzleMixer(To6, Po6, beta, f, fab, Po2, To2)
Prnm = 0.8;
To7 = (beta*To2+(1+f+fab)*To6)/(1+beta+f+fab);
gamma = 1.44 - (1.39*10^-4)*To7 + (3.57*10^-8)*To7;
Po7 = Po6*Prnm*((Po2/Po6)^(beta/(1+f+fab)))*((To7/To6)^(gamma/(gamma-1)))*((To6/To2)^((gamma*beta)/(gamma-1)/(1+f+fab)));
end
<<<<<<< HEAD

function [Tec, uec] = combinedNozzle(To7, Po7, Pa)
=======
%+
function [Tec, uec] = combinedNozzle(To7, Po7, Pa, MW)
>>>>>>> c2925fbee04d1c1f09affd74e3c9d272a73fc902
gamma = 1.37;
nc = 0.95;
Cp = gamma*(R_/MW)/(gamma-1);
Tec = To7*(1-nc*(1-(Pa/Po7)^((gamma-1)/gamma)));
uec = sqrt(2*Cp*(To7-Tec));
<<<<<<< HEAD
end
=======
end
%%

>>>>>>> c2925fbee04d1c1f09affd74e3c9d272a73fc902

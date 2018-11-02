%% AE 4451 Jet & Rocket Propulsion Project
%10/26/2018
%Authors: Loren Isakson, Brandon Campanile, Matthew Yates
clear;clc
close all
%% Inputs

%NOTE: Always use Pa, K, kg, m, sec, J, N, etc for consistency or I will
%find you - Loren

%flight conditions
Ta = 220; %ambient temperature
Pa = 10*10^3; %ambient pressure
M = 1.5; %flight mach number
u = M*sqrt(1.4*8314*Ta/28.8);

%pressure ratios
Prf = 1.2; %fan
Prc = 30; %compressor
Prb = 0.98; %burner
Prab = 0.97; %afterburner

%fuel/air mass ratios
f = 0.018; %main burner fuel-air ratio
fab = 0.01; %afterburner fuel-air ratio
beta = 2; %bypass ratio
b = 0.1; %bleed ratio

%specific heat ratios
yd = 1.4; %diffuser
yf = 1.4; %fan
yc = 1.38; %compressor
yb = 1.33; %burner
yt = 1.33; %turbine
ytm = 1.34; %turbine mixer
yft = 1.33; %fan turbine
yab = 1.32; %afterburner
yn = 1.35; %core nozzle
yfn = 1.4; %fan nozzle
ycn = 1.37; %combined nozzle

%molecular weights
MWd = 28.8; %diffuser
MWf = 28.8; %fan
MWc = 28.8; %compressor
MWb = 28.8; %burner
MWt = 28.8; %turbine
MWtm = 28.8; %turbine mixer
MWft = 28.8; %fan turbine
MWab = 28.8; %afterburner
MWn = 28.8; %core nozzle
MWfn = 28.8; %fan nozzle
MWnm = 28.8; %nozzle mixer
MWcn = 28.8; %combined nozzle

%efficiencies
nd = 0.92; %diffuser
nf = 0.90; %fan
nc = 0.90; %compressor
nb = 0.99; %burner
nt = 0.92; %turbine
nft = 0.92; %fan turbine
nab = 0.96; %afterburner
nn = 0.95; %core nozzle
nfn = 0.97; %fan nozzle
ncn = 0.95; %combined nozzle
nfp = 0.35; %fuel pump

%% test case

[To1, Po1] = diffuser(Ta, Pa, M, yd, nd);
[To2, Po2, wf_ma] = fan(To1, Po1, Prf, beta, MWf, yf, nf);
[To3, Po3, wc_ma] = compressor(To2, Po2, Prc, Prf, MWc, yc, nc);
[To4, Po4, wp_ma, fmax] = burner(To3, Po3, Prb, b, f, MWb, yb, nb, nfp);
[To5_1, Po5_1] = turbine(To4, Po4, wc_ma, wp_ma, b, f, MWt, yt, nt);
[To5_m, Po5_m] = turbineMixer(To5_1, Po5_1, To3, f, b, MWtm, ytm);
[To5_2, Po5_2] = fanTurbine(To5_m, Po5_m, wf_ma, f, MWft, yft, nft);
[To6, Po6, fabmax] = afterburner(Po5_2, To5_2, Prab, f, fab, fmax, MWab, yab, nab);
[To7, Po7, ynm] = nozzleMixer(To6, Po6, beta, f, fab, Po2, To2);
[Tec, Pec, uec] = combinedNozzle(To7, Po7, Pa, MWcn, ycn, ncn);

[Tef, Pef, uef] = fanNozzle(To2, Po2, Pa, MWfn, yfn, nfn);
[Te, Pe, ue] = coreNozzle(To6, Po6, Pa, MWcn, yn, nn);

[ST, TSFC, effth, effp, effo] = performance(f, fab, uec, uec, u, beta, Pa, M);
[STf, TSFCf, effthf, effpf, effof] = performance(f, fab, ue, uef, u, beta, Pa, M);
titletop = [{'Inputs'}, cell(1,11)];
inputs = {'Ta(K)','Pa(kPa)','M','Prc','Prf','beta','b','f','fab';
           Ta,Pa/1000,M,Prc,Prf,beta,b,f,fab};
inputs = [inputs cell(2,3)];
titlebot = [{'Ouputs'}, cell(1,11)];
outputs = {'To1(K)','Po1(kPa)','To2(K)','Po2(kPa)','To3(K)','Po3(kPa)','To4(K)','Po4(kPa)','To5.1','Po5.1(kPa)','To5.m(K)','Po5.m(kPa)';
            To1,Po1/1000,To2,Po2/1000,To3,Po3/1000,To4,Po4/1000,To5_1,Po5_1/1000,To5_m,Po5_m/1000;
           'To5.2(K)','Po5.2(kPa)','To6(K)','Po6(kPa)','Te(K)','Pe(K)','Tef(K)','Pef(kPa)','To7(K)','ynm','Po7(kPa)','Tec(K)';
            To5_2,Po5_2/1000,To6,Po6/1000,Te,Pe/1000,Tef,Pef,To7,ynm,Po7,Tec};

perform = {'ue(m/s)','uef(m/s)','ST(kNs/kg)','TSFC(kg/kNs)','nth(%)','no(%)','Wc(kJ/kg)','Wp(kJ/kg)','Wft(kJ/kg)','fmax','fmaxab';
            ue,uef,ST/1000,TSFC,effth*100,effo*100,wc_ma/1000,wp_ma/1000,wf_ma/1000,fmax,fabmax};
perform = [perform, cell(2,1)];
warning('off','MATLAB:xlswrite:AddSheet')
xlswrite('Test_Case.xlsx',[titletop;inputs;cell(1,12);titlebot;outputs;perform]);
disp('Run Completed');

%% Optimization

function table = engineAnalysis(Ta, Pa, Pf, M, Prc, Prf, beta, b, f, fab, ab, Nmix, Tmax_ab, Prab)
% Ta = ambient temperature [Kelvin]; Pa = ambient pressure [kPa]; Pf = fuel storage pressure [kPa]; M = flight mach number
% Prc = compressor stagnation pressure ratio; Prf = fan stagnation; pressure; ratio; B = bypass ratio; b = bleed ratio
% f = main burner fuel-air ratio; fab = afterburner fuel-air ratio; ab = afterburner boolean; Nmix = nozzle mixing boolean
% Tmax_ab = max stagnation temp afternburner; Prab = stagnation pressure ratio afterburner
syms Te Tef To7 Tec Po7 uec uef ue Po6 To6
MW=[MWf, MWc, MWb, MWt1, MWtm, MWt2, MWcN, MWfn, MWn, MWab];
u=M*sqrt(1.4*Ta*R_/MW(0));
ST=0;
i=0;
while 2.73>ST || ST>2.75 % Loop until ST=2.74 to meet goal 1. Then keeping this value for Prf, vary other parameters to maximize ST. 
    % while loop is temporary until Prf and engine cycle is determined for each case. then change to an optimization thing.
    if i <=.4
        Prf=1.1+i;
    else
        error('get real')
    end
    i=i+.005;
    Prc=60/Prf;
    % everything above to while loop is temporary
    if ~Prf % add user input logic
        % Turbofan
        [To1, Po1] = diffuser(Ta, Pa, M);
        [To2, Po2, wf_ma] = fan(To1, Po1, Prf, beta, MW(0));
        [To3, Po3, wc_ma] = compressor(To2, Po2, Prc, Prf, MW(1));
        [To4, Po4, wp_ma] = burner(Po3, b, f, MW(2));
        [To5_1, Po5_1] = turbine(To4, Po4, wc_ma, wp_ma, b, f, MW(3));
        [To5_m, Po5_m] = turbineMixer(To5_1, Po5_1, To3, f, b, MW(4));
        [To5_2, Po5_2] = fanTurbine(To5_m, Po5_m, wf_ma, f, MW(5));
        if ab && Nmix
            [To6, Po6] = afterburner(Po5_2, To5_2, Prab, f, fab, MW(9));
            [To7, Po7] = nozzleMixer(To6, Po6, beta, f, fab, Po2, To2);
            [Tec, uec] = combinedNozzle(To7, Po7, Pa, MW(6));
            [ST, TSFC, effth, effp, effo] = performance(f, fab, uec, uec, u, beta, Pa, M);
        elseif ab && ~Nmix
            [To6, Po6] = afterburner(Po5_2, To5_2, Prab, f, fab, MW(9));
            [Tef, uef] = fanNozzle(To2, Po2, Pa, MW(7));
            [Te, ue] = coreNozzle(To6, Po6, Pa, MW(8));
            [ST, TSFC, effth, effp, effo] = performance(f, fab, ue, uef, u, beta, Pa, M);
        elseif ~ab && Nmix
            [To7, Po7] = nozzleMixer(To5_2, Po5_2, beta, f, fab, Po2, To2);
            [Tec, uec] = combinedNozzle(To7, Po7, Pa, MW(6));
            [ST, TSFC, effth, effp, effo] = performance(f, 0, uec, uec, u, beta, Pa, M);
        else
            [Te, ue] = coreNozzle(To5_2, Po5_2, Pa, MW(9));
            [ST, TSFC, effth, effp, effo] = performance(f, 0, ue, ue, u, beta, Pa, M);
        end
    else % turbojet
        [To1, Po1] = diffuser(Ta, Pa, M);
        [To3, Po3, wc_ma] = compressor(To1, Po1, Prc, Prf, MW(1));
        [To4, Po4, wp_ma] = burner(Po3, b, f, MW(2));
        [To5_1, Po5_1] = turbine(To4, Po4, wc_ma, wp_ma, b, f, MW(3));
        [To5_m, Po5_m] = turbineMixer(To5_1, Po5_1, To3, f, b, MW(4));
        if ab
            [To6, Po6] = afterburner(Po5_2, To5_2, Prab, f, fab, MW(9));
            [Te, ue] = coreNozzle(To6, Po6, Pa, MW(8));
            [ST, TSFC, effth, effp, effo] = performance(f, fab, ue, ue, u, 0, Pa, M);
        else
            [Te, ue] = coreNozzle(To5_m, Po5_m, Pa, MW(8));
            [ST, TSFC, effth, effp, effo] = performance(f, 0, ue, ue, u, 0, Pa, M);
        end
    end
    warning('off','MATLAB:xlswrite:AddSheet')
    table = xlswrite('Results.xls', ['To1', 'To2', 'To3', 'To4', 'To5_1', 'To5_m', 'To5_2', 'To6', 'To7';
        To1, To2, To3, To4, To5_1, To5_m, To5_2, To6, To7;
        'Po1', 'Po2', 'Po3', 'Po4', 'Po5_1', 'Po5_m', 'Po5_2', 'Po6', 'Po7';
        Po1, Po2, Po3, Po4, Po5_1, Po5_m, Po5_2, Po6, Po7;
        'Te_combined', 'ue_combined', 'Te_fan', 'ue_fan', 'Te_core', 'ue_core';
        Tec, uec, Tef, uef, Te, ue;
        'ST', 'TSFC', 'eff_th', 'eff_p', 'eff_o';
        ST, TSFC, effth, effp, effo]);
end
end

%% Component Functions

function [To1, Po1] = diffuser(Ta, Pa, M, yd, nd)
%static ambient temp, press, mach number, gamma, adiabatic efficiency
gamma = yd;
if M<=1
    rd=1;
else
    rd=1-.075*(M-1)^1.35;
end
To1 = Ta.*(1+0.5.*(gamma-1).*M.^2);
Po1 = (Pa.*(1+nd.*(To1/Ta - 1)).^(gamma/(gamma-1)))*rd;
end

function [To2, Po2, wf_ma] = fan(To1, Po1, Prf, beta, MW, yf, nf)
%stagnation temp, press, pressure ratio, polytropic efficiency
R = 8314;
gamma = yf;
Cp = gamma*(R/MW)/(gamma-1);
if Prf < 1.1 || Prf > 1.5
    error('Fan pressure ratio is constrained 1.1 <= Prf <= 1.5');
end 
To2 = To1.*(Prf).^((gamma-1)/gamma/nf);
Po2 = Po1.*Prf;
wf_ma = Cp*(1+beta)*(To2-To1); %work done by the fan on the fluid
end

function [To3, Po3, wc_ma] = compressor(To2, Po2, Prc, Prf, MW, yc, nc)
gamma = yc;
R = 8314;
Cp = gamma*(R/MW)/(gamma-1);
if Prc > 60/Prf
    error('Compressor pressure ratio must be less than 60/fan pressure ratio');
end
Po3 = Po2*Prc;
To3 = To2*(Prc).^((gamma-1)/gamma/nc);
wc_ma = Cp*(To3-To2);
end

function [To4, Po4, wp_ma, fmax] = burner(To3, Po3, Prb, b, f, MW, yb, nb, nfp)
gamma = yb;
R = 8314;
HVf = 45e6;
Cp = gamma*(R/MW)/(gamma-1);
Tomax = 1300; %kelvin
Cb = 700; %kelvin
bmax = 0.12;
Tmax = Tomax + Cb*(b/bmax)^0.5;
fmax = (1-b)*(1-To3/Tmax)/((nb*HVf/Cp/Tmax)-1);
if f > fmax
    error('f is greater than fmax = %d',fmax);
end
To4 = (1/(1-b+f))*((1-b)*To3+nb*HVf*f/Cp);
Po4 = Po3*Prb;
%Fuel Pump
Pf = 104*10^3; %fuel storage pressure
deltaP_inject = 550*10^3;
np = nfp;
f_density = 780;
wp_ma = f*deltaP_inject/np/f_density;
end

function [To5_1, Po5_1] = turbine(To4, Po4, wc_ma, wp_ma, b, f, MW, yt, nt)
npt = nt;
gamma = yt;
R = 8314;
Cp = gamma*(R/MW)/(gamma-1);
To5_1 = To4 - (wc_ma + wp_ma)/(Cp*(1+f-b));
TR = To5_1/To4;
Po5_1 = Po4*(TR^(1/npt))^(gamma/(gamma-1));
end

function [To5_m, Po5_m] = turbineMixer(To5_1, Po5_1, To3, f, b, MW, ytm)
gamma = ytm;
R = 8314;
Cp = gamma*(R/MW)/(gamma-1);
To5_m = (b*To3 + (1+f-b)*To5_1)/(1+f);
Po5_m = Po5_1*((To5_m/To5_1)^(gamma/(gamma-1)))*(To5_1/To3)^(gamma*b/((gamma-1)*(1+f)));
end

function [To5_2, Po5_2] = fanTurbine(To5_m, Po5_m, wf_ma, f, MW, yft, nft)
gamma = yft;
R = 8314;
npft = nft;
Cp = gamma*(R/MW)/(gamma-1);
To5_2 = To5_m - wf_ma/(Cp*(1+f));
TR = To5_2/To5_m;
Po5_2 = Po5_m*(TR^(1/npft))^(gamma/(gamma-1));
end

function [To6, Po6, fabmax] = afterburner(Po5_2, To5_2, Prab, f, fab, fmax, MW, yab, nab)
gamma = yab;
R = 8314;
HVf = 45e6;
cpab = gamma*(R/MW)/(gamma-1);
Tmax_ab = 2200;
fabmax = (1+fmax)*((Tmax_ab/To5_2)-1)/((nab*HVf/cpab/To5_2)-(Tmax_ab/To5_2));
if fab > fabmax
    error('fab is greater than fabmax = %d',fabmax);
end
if fab>0
    PR = Prab;
else
    PR = 1;
end
Po6 = Po5_2*PR;
%To6 = min((1/(1+f+fab))*((1+f)*To5_2+nab*HVf*fab/cpab), Tmax_ab)
To6 = (1/(1+f+fab))*((1+f)*To5_2+nab*HVf*fab/cpab);
end

function [Te, Pe, ue] = coreNozzle(To6, Po6, Pa, MW, yn, nn)
Pe = Pa;
gamma = yn;
R = 8314;
nc = nn;
Cp = gamma*(R/MW)/(gamma-1);
Te = To6*(1-nc*(1-(Pe/Po6)^((gamma-1)/gamma)));
ue = sqrt(2*Cp*(To6-Te));
end

function [Tef, Pef, uef] = fanNozzle(To2, Po2, Pa, MW, yfn, nfn)
Pef = Pa;
gamma = yfn;
R = 8314;
nf = nfn;
Cp = gamma*(R/MW)/(gamma-1);
Tef = To2*(1-nf*(1-(Pa/Po2)^((gamma-1)/gamma)));
uef = sqrt(2*Cp*(To2-Tef));
end

function [To7, Po7, ynm] = nozzleMixer(To6, Po6, beta, f, fab, Po2, To2)
Prnm = 0.8;
To7 = (beta*To2+(1+f+fab)*To6)/(1+beta+f+fab);
gamma = 1.44 - (1.39*10^-4)*To7 + (3.57*10^-8)*To7^2;
ynm = gamma;
Po7 = Po6*Prnm*((Po2/Po6)^(beta/(1+beta+f+fab)))*((To7/To6)^(gamma/(gamma-1)))*((To6/To2)^((gamma*beta)/(gamma-1)/(1+beta+f+fab)));
end

function [Tec, Pec, uec] = combinedNozzle(To7, Po7, Pa, MW, ycn, ncn)
Pec = Pa;
gamma = ycn;
R = 8314;
nc = ncn;
Cp = gamma*(R/MW)/(gamma-1);
Tec = To7*(1-nc*(1-(Pa/Po7)^((gamma-1)/gamma)));
uec = sqrt(2*Cp*(To7-Tec));
end

%% Specific Thrust, Specific Fuel Consumption, Efficiencies

function [ST, TSFC, effth, effp, effo] = performance(f, fab, ue, uef, u, beta, Pa, M)
Patm = 101.3*10^3;
CB = 245;
deltaD = CB*M^2*Pa/Patm*beta^1.5;
HVf = 45e6;
ST = (1+f+fab)*ue+beta*uef-(beta+1)*u - deltaD;
TSFC = (f+fab)/ST;
effth = (beta*uef^2+(1+f+fab)*ue^2-(beta+1)*u^2)/(2*(f+fab)*HVf);
effp = 2*ST*u/(beta*uef^2+(1+f+fab)*ue^2-(beta+1)*u^2);
effo = effth*effp;
end

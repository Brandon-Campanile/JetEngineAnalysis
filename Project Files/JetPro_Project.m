%% AE 4451 Jet & Rocket Propulsion Project
%10/26/2018
%Authors: Loren Isakson, Brandon Campanile, Matthew Yates


%% Main Function

function out = JetPro_Project(T, eType, Nmix, Ta, Pa, Pf, M, Prf, Prc, Prb, Prab, Prnm, beta, b, f, fab, Tomax, Tmax_ab, MW, eff, y, HVf)
% Ta = ambient temperature [Kelvin]; Pa = ambient pressure [kPa]; Pf = fuel storage pressure [kPa];
% M = flight mach number Prc = compressor stagnation pressure ratio; Prf = fan stagnation pressure ratio;
% Prb = burner stagnation pressure ratio; beta = bypass ratio; b = bleed ratio; f = main burner fuel-air ratio;
% fab = afterburner fuel-air ratio; Nmix = nozzle mixing boolean; Tmax = burner max stagnation temp [K];
% Tmax_ab = max stagnation temp afternburner [K]; Prab = stagnation pressure ratio afterburner; MW = list of all
% molecular weights; y = list of all specific heat ratios; eff = list of component efficiencies
% eType = 1 for turbofan 0 for turbojet; Nmix = nozzle mixing?; T = 1 for final run, 0 for optimization run 

ynm=y(8);
R = 8314;
u = M*sqrt(y(1)*Ta*R/MW(1)); %vaircraft velocity

% Run components based on engine cycles
if strcmp(eType,'Turbojet')
    beta=0;
    [To1, Po1] = diffuser(Ta, Pa, M, y(1), eff(1));
    [To3, Po3, wc_ma] = compressor(To1, Po1, Prc, MW(3), y(3), eff(3));
    [To4, Po4, wp_ma, fmax, f] = burner(To3, Po3, Prb, b, f, MW(4), y(4), eff(4), eff(11), HVf, Tomax, Pf, fab);
    [To5_1, Po5_1] = turbine(To4, Po4, wc_ma, wp_ma, b, f, MW(5), y(5), eff(5));
    [To5_m, Po5_m] = turbineMixer(To5_1, Po5_1, To3, f, b, MW(6), y(6));    
    Po1=Po1/1000;
    Po2='NA';
    Po3=Po3/1000;
    wf_ma='NA';
    To2='NA';
    Po5_2='NA';
    To5_2='NA';
    Po7='NA';
    To7='NA';
    Tef='NA';
    Pef='NA';
    uef='NA';
    Tec='NA';
    Pec='NA';
    uec='NA';
    Po4=Po4/1000;
    wc_ma=wc_ma/1000;
    wp_ma=wp_ma/1000;
    Po5_1=Po5_1/1000;
    if fab>0 % with afterburner
        [To6, Po6, fabmax, fab] = afterburner(Po5_m, To5_m, Prab, f, fab, fmax, MW(8), y(8), eff(7), HVf, Tmax_ab);
        [Te, Pe, ue] = coreNozzle(To6, Po6, Pa, MW(9), y(9), eff(8));
        [ST, TSFC, effth, effp, effo] = performance(f, fab, ue, 0, u, beta, Pa, M, HVf);
        Po5_m=Po5_m/1000;
        Po6=Po6/1000;
        Pe=Pe/1000;
    else % without afterburner
        [Te, Pe, ue] = coreNozzle(To5_m, Po5_m, Pa, MW(9), y(9), eff(8));
        [ST, TSFC, effth, effp, effo] = performance(f, fab, ue, 0, u, beta, Pa, M, HVf);
        Po5_m=Po5_m/1000;
        Po6='NA';
        To6='NA';
        Pe=Pe/1000;
        fabmax='NA';
    end
elseif strcmp(eType,'Turbofan')
    [To1, Po1] = diffuser(Ta, Pa, M, y(1), eff(1));
    [To2, Po2, wf_ma] = fan(To1, Po1, Prf, beta, MW(2), y(2), eff(2));
    [To3, Po3, wc_ma] = compressor(To2, Po2, Prc, MW(3), y(3), eff(3));
    [To4, Po4, wp_ma, fmax, f] = burner(To3, Po3, Prb, b, f, MW(4), y(4), eff(4), eff(11), HVf, Tomax, Pf, fab);
    [To5_1, Po5_1] = turbine(To4, Po4, wc_ma, wp_ma, b, f, MW(5), y(5), eff(5));
    [To5_m, Po5_m] = turbineMixer(To5_1, Po5_1, To3, f, b, MW(6), y(6));
    [To5_2, Po5_2] = fanTurbine(To5_m, Po5_m, wf_ma, f, MW(7), y(7), eff(6));
    Po1=Po1/1000;
    Po3=Po3/1000;
    Po4=Po4/1000;
    wf_ma=wf_ma/1000;
    wc_ma=wc_ma/1000;
    wp_ma=wp_ma/1000;
    Po5_1=Po5_1/1000;
    Po5_m=Po5_m/1000;
    if fab>0 && Nmix % with afterburner and nozzle mixing 
        ue='NA';
        Tef='NA';
        Pef='NA';
        uef='NA';
        Te='NA';
        Pe='NA';
        ue='NA';
        [To6, Po6, fabmax, fab] = afterburner(Po5_2, To5_2, Prab, f, fab, fmax, MW(8), y(8), eff(7), HVf, Tmax_ab);
        [To7, Po7, ynm] = nozzleMixer(To6, Po6, beta, f, fab, Po2, To2, Prnm);
        [Tec, Pec, uec] = combinedNozzle(To7, Po7, Pa, MW(12), y(11), eff(10));
        [ST, TSFC, effth, effp, effo] = performance(f, fab, uec, uec, u, beta, Pa, M, HVf);
        Po7=Po7/1000;
        Po6=Po6/1000;
        Po5_2=Po5_2/1000;
        Pec=Pec/1000;
        Po2=Po2/1000;
    elseif fab>0 && ~Nmix % with afterburner, no nozzle mixing
        To7='NA';
        Po7='NA';
        ynm='NA';
        uec='NA';
        Tec='NA';
        Pec='NA';
        [To6, Po6, fabmax, fab] = afterburner(Po5_2, To5_2, Prab, f, fab, fmax, MW(8), y(8), eff(7), HVf, Tmax_ab);
        [Tef, Pef, uef] = fanNozzle(To2, Po2, Pa, MW(10), y(10), eff(9));
        [Te, Pe, ue] = coreNozzle(To6, Po6, Pa, MW(9), y(9), eff(8));
        [ST, TSFC, effth, effp, effo] = performance(f, fab, ue, uef, u, beta, Pa, M, HVf);
        Pef = Pef/1000;
        Pe=Pe/1000;
        Po5_2=Po5_2/1000;
        Po6=Po6/1000;
        Po2=Po2/1000;
    elseif fab==0 && Nmix % with nozzle mixing, no afterburner
        To6='NA';
        Po6='NA';
        ue='NA';
        uef='NA';
        Pe='NA';
        Te='NA';
        Tef='NA';
        Pef='NA';
        fabmax='NA';
        [To7, Po7, ynm] = nozzleMixer(To5_2, Po5_2, beta, f, fab, Po2, To2, Prnm);
        [Tec, Pec, uec] = combinedNozzle(To7, Po7, Pa, MW(12), y(11), eff(10));
        [ST, TSFC, effth, effp, effo] = performance(f, fab, uec, uec, u, beta, Pa, M, HVf);
        Po7=Po7/1000;
        Pec=Pec/1000;
        Po5_2=Po5_2/1000;
        Po7=Po7/1000;
        Po2=Po2/1000;
    else % no nozzle mixing or afterburner
        To6='NA';
        Po6='NA';
        To7='NA';
        Po7='NA';
        uec='NA';
        Tec='NA';
        Pec='NA';
        ynm='NA';
        fabmax='NA';
        [Tef, Pef, uef] = fanNozzle(To2, Po2, Pa, MW(10), y(10), eff(9));
        [Te, Pe, ue] = coreNozzle(To5_2, Po5_2, Pa, MW(9), y(9), eff(8));
        [ST, TSFC, effth, effp, effo] = performance(f, fab, ue, uef, u, beta, Pa, M, HVf);
        Pef = Pef/1000;
        Pe = Pe/1000;
        Po5_2 = Po5_2/1000;
        Po2=Po2/1000;
        Pef=Pef/1000;
        Pe=Pe/1000;
    end
else %Ramjet
    beta=0;
    [To1, Po1] = diffuser(Ta, Pa, M, y(1), eff(1));
    To2='NA';
    Po2='NA';
    To3='NA';
    Po3='NA';
    wc_ma='NA';
    [To4, Po4, wp_ma, fmax, f] = burner(To1, Po1, Prb, b, f, MW(4), y(4), eff(4), eff(11), HVf, Tomax, Pf, fab);
    wp_ma=wp_ma/1000;
    To5_1='NA';
    Po5_1='NA';
    To5_m='NA';
    Po5_m='NA';
    To5_2='NA';
    Po5_2='NA';
    wf_ma='NA';
    Tec='NA';
    Pec='NA';
    uec='NA';
    Po1=Po1/1000;
    if fab>0 % with afterburner
        To7 = 'NA';
        Po7='NA';
        Tef='NA';
        Pef='NA';
        uef='NA';
        [To6, Po6, fabmax, fab] = afterburner(Po4, To4, Prab, f, fab, fmax, MW(8), y(8), eff(7), HVf, Tmax_ab);
        [Te, Pe, ue] = coreNozzle(To6, Po6, Pa, MW(9), y(9), eff(8));
        [ST, TSFC, effth, effp, effo] = performance(f, fab, ue, 0, u, beta, Pa, M, HVf);
        Po4=Po4/1000;
        Po6=Po6/1000;
        Pe=Pe/1000;
    else % without afterburner
        To6='NA';
        Po6='NA';
        To7='NA';
        Po7='NA';
        Tef='NA';
        Pef='NA';
        uef='NA';
        fabmax='NA';
        [Te, Pe, ue] = coreNozzle(To4, Po4, Pa, MW(9), y(9), eff(8));
        [ST, TSFC, effth, effp, effo] = performance(f, fab, ue, 0, u, beta, Pa, M, HVf);
        Po4=Po4/1000;
        Pe=Pe/1000;
    end
end
if T
    % Generate Table for outputs
    titletop = [{'Inputs'}, cell(1,11)];
    inputs = {'Ta(K)','Pa(kPa)','M','Prc','Prf','beta','b','f','fab';
        Ta,Pa/1000,M,Prc,Prf,beta,b,f,fab};
    inputs = [inputs cell(2,3)];
    titlebot = [{'Outputs'}, cell(1,11)];
    output1 = {'To1(K)','Po1(kPa)','To2(K)','Po2(kPa)','To3(K)','Po3(kPa)','To4(K)','Po4(kPa)','To5.1','Po5.1(kPa)','To5.m(K)','Po5.m(kPa)';
        To1,Po1,To2,Po2,To3,Po3,To4,Po4,To5_1,Po5_1,To5_m,Po5_m};
    output2 = {'To5.2(K)','Po5.2(kPa)','To6(K)','Po6(kPa)','To7(K)','Po7(kPa)','Te(K)','Pe(kPa)','ue(m/s)','Tef(K)','Pef(kPa)','uef(m/s)';
        To5_2,Po5_2,To6,Po6,To7,Po7,Te,Pe,ue,Tef,Pef,uef};
    
    perform1 = {'Tec(K)','Pec(kPa)','uec(m/s)','ynm','ST(kNs/kg)','TSFC(kg/kNs)','nth(%)','np(%)','no(%)','Wc(kJ/kg)','Wp(kJ/kg)','Wft(kJ/kg)';
        Tec, Pec, uec, ynm,ST/1000,TSFC*1000,effth*100,effp*100,effo*100,wc_ma,wp_ma,wf_ma};
    perform2 = {'fmax','fmaxab','','','','','','','','','','';fmax,fabmax,'','','','','','','','','',''};
    warning('off','MATLAB:xlswrite:AddSheet')
    xlswrite('Results.xlsx',[titletop;inputs;cell(2,12);titlebot;output1;cell(1,12);output2;cell(1,12);perform1;cell(1,12);perform2]);
end

if strcmp(eType,'Turbojet')
    out = [ST, TSFC, To3, To5_m, f, fab];
elseif strcmp(eType,'Turbofan')
    out = [ST, TSFC, To3, To5_2, f, fab];
else
    out = [ST, TSFC, To1, To4, f, fab];
end
end

%% Component Functions

function [To1, Po1] = diffuser(Ta, Pa, M, yd, nd)
%static ambient temp, press, mach number, gamma, adiabatic efficiency

if M<=1
    rd=1;
else
    rd=1-.075*(M-1)^1.35;
end
To1 = Ta.*(1+0.5.*(yd-1).*M.^2);
Po1 = (Pa.*(1+nd.*(To1/Ta - 1)).^(yd/(yd-1)))*rd;
end

function [To2, Po2, wf_ma] = fan(To1, Po1, Prf, beta, MW, yf, nf)
%stagnation temp, press, pressure ratio, polytropic efficiency
R = 8314;
Cp = yf*(R/MW)/(yf-1);
To2 = To1.*(Prf).^((yf-1)/yf/nf);
Po2 = Po1.*Prf;
wf_ma = Cp*(1+beta)*(To2-To1); % work done by the fan on the fluid
end

function [To3, Po3, wc_ma] = compressor(To2, Po2, Prc, MW, yc, nc)
R = 8314;
Cp = yc*(R/MW)/(yc-1);

Po3 = Po2*Prc;
To3 = To2*(Prc).^((yc-1)/yc/nc);
wc_ma = Cp*(To3-To2);
end

function [To4, Po4, wp_ma, fmax, f] = burner(To3, Po3, Prb, b, f, MW, yb, nb, nfp, HVf, Tomax, Pf, fab)
Cb = 700; %kelvin
R = 8314;
bmax = 0.12;
deltaP_inject = 550*10^3;
f_density = 780;

Tmax = Tomax + Cb*(b/bmax)^0.5;
Cp = yb*(R/MW)/(yb-1);
fmax = (1-b)*(1-To3/Tmax)/((nb*HVf/Cp/Tmax)-1);

To4 = (1/(1-b+f))*((1-b)*To3+nb*HVf*f/Cp);
Po4 = Po3*Prb;

%Fuel Pump
wp_ma = (f+fab)*(Po3+deltaP_inject-Pf)/nfp/f_density;
end

function [To5_1, Po5_1] = turbine(To4, Po4, wc_ma, wp_ma, b, f, MW, yt, nt)
R = 8314;
Cp = yt*(R/MW)/(yt-1);
To5_1 = To4 - (wc_ma + wp_ma)/(Cp*(1+f-b));
TR = To5_1/To4;
Po5_1 = Po4*(TR^(1/nt))^(yt/(yt-1));
end

function [To5_m, Po5_m] = turbineMixer(To5_1, Po5_1, To3, f, b, MW, ytm)
%R = 8314;
%Cp = ytm*(R/MW)/(ytm-1);
To5_m = (b*To3 + (1+f-b)*To5_1)/(1+f);
Po5_m = Po5_1*((To5_m/To5_1)^(ytm/(ytm-1)))*(To5_1/To3)^(ytm*b/((ytm-1)*(1+f)));
end

function [To5_2, Po5_2] = fanTurbine(To5_m, Po5_m, wf_ma, f, MW, yft, nft)
R = 8314;
Cp = yft*(R/MW)/(yft-1);
To5_2 = To5_m - wf_ma/(Cp*(1+f));
TR = To5_2/To5_m;
Po5_2 = Po5_m*(TR^(1/nft))^(yft/(yft-1));
end

function [To6, Po6, fabmax, fab] = afterburner(Po5_2, To5_2, Prab, f, fab, fmax, MW, yab, nab, HVf, Tmax_ab)
R = 8314;
cpab = yab*(R/MW)/(yab-1);
fabmax = (1+fmax)*((Tmax_ab/To5_2)-1)/((nab*HVf/cpab/To5_2)-(Tmax_ab/To5_2));

Po6 = Po5_2*Prab;
To6 = (1/(1+f+fab))*((1+f)*To5_2+nab*HVf*fab/cpab);
end

function [Te, Pa, ue] = coreNozzle(To6, Po6, Pa, MW, yn, nn)
R = 8314;
Cp = yn*(R/MW)/(yn-1);
Te = To6*(1-nn*(1-(Pa/Po6)^((yn-1)/yn)));
ue = sqrt(2*Cp*(To6-Te));
if ~isreal(ue)
    ue=0;
end
end

function [Tef, Pa, uef] = fanNozzle(To2, Po2, Pa, MW, yfn, nfn)
R = 8314;
Cp = yfn*(R/MW)/(yfn-1);
Tef = To2*(1-nfn*(1-(Pa/Po2)^((yfn-1)/yfn)));
uef = sqrt(2*Cp*(To2-Tef));
end

function [To7, Po7, gamma] = nozzleMixer(To6, Po6, beta, f, fab, Po2, To2, Prnm)
To7 = (beta*To2+(1+f+fab)*To6)/(1+beta+f+fab);
gamma = 1.44 - (1.39*10^-4)*To7 + (3.57*10^-8)*To7^2;
Po7 = Po6*Prnm*((Po2/Po6)^(beta/(1+beta+f+fab)))*((To7/To6)^(gamma/(gamma-1)))*((To6/To2)^((gamma*beta)/(gamma-1)/(1+beta+f+fab)));
end

function [Tec, Pa, uec] = combinedNozzle(To7, Po7, Pa, MW, ycn, ncn)
R = 8314;
Cp = ycn*(R/MW)/(ycn-1);
Tec = To7*(1-ncn*(1-(Pa/Po7)^((ycn-1)/ycn)));
uec = sqrt(2*Cp*(To7-Tec));
if ~isreal(uec)
    uec=0;
end
end


%% Specific Thrust, Specific Fuel Consumption, Efficiencies

function [ST, TSFC, effth, effp, effo] = performance(f, fab, ue, uef, u, beta, Pa, M, HVf)
Patm = 101.3*10^3;
CB = 245;
deltaD = CB*M^2*Pa/Patm*beta^1.5;
ST = (1+f+fab)*ue+beta*uef-(beta+1)*u - deltaD;
TSFC = (f+fab)/ST;
effth = (beta*uef^2+(1+f+fab)*ue^2-(beta+1)*u^2)/(2*(f+fab)*HVf);
effp = 2*ST*u/(beta*uef^2+(1+f+fab)*ue^2-(beta+1)*u^2);
effo = effth*effp;
end

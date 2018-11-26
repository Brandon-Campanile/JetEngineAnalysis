function ST_max = Maximize(ab, eType, Nmix, Ta, Pa, Pf, M, Prf, Prc, Prb, Prab, Prnm, beta, ~, ~, ~, Tomax, Tmax_ab, MW, eff, y, HVf)

T = 0;

if ab % w/ afterburner
        % [fuel-air ratio, afterburner fuel-air ratio, bleed ratio]
        lb = [.001, .0005, 0]; % lower bound
        ub = [.1, .1, .12];        % upper bound
        x0 = [.05, .05, .05];
else
        lb = [.001, 0, 0];
        ub = [.1, 0, .12];
        x0 = [.05, 0, .05];
end

func = @(x)singleOut(x, T, eType, Nmix, Ta, Pa, Pf, M, Prf, Prc, Prb, Prab, Prnm, beta, Tomax, Tmax_ab, MW, eff, y, HVf);

nlc = @(x)nonlcon(x, T, eType, Nmix, Ta, Pa, Pf, M, Prf, Prc, Prb, Prab, Prnm, beta, Tomax, Tmax_ab, MW, eff, y, HVf);

problem = createOptimProblem('fmincon','objective', func,'x0',x0,'lb',lb,'ub',ub,'nonlcon', nlc);

gs=GlobalSearch;

ST_max = run(gs,problem);

end

function ST = singleOut(x, T, eType, Nmix, Ta, Pa, Pf, M, Prf, Prc, Prb, Prab, Prnm, beta, Tomax, Tmax_ab, MW, eff, y, HVf)
f = x(1);
fab = x(2);
b = x(3);
output = JetPro_Project(T, eType, Nmix, Ta, Pa, Pf, M, Prf, Prc, Prb, Prab, Prnm, beta, b, f, fab, Tomax, Tmax_ab, MW, eff, y, HVf);
ST = -output{1};
end

function [g,ceq] = nonlcon(x, T, eType, Nmix, Ta, Pa, Pf, M, Prf, Prc, Prb, Prab, Prnm, beta, Tomax, Tmax_ab, MW, eff, y, HVf)
R=8314;
CB=700;
bmax=.12;
f = x(1);
fab = x(2);
b = x(3);

out2 = JetPro_Project(T, eType, Nmix, Ta, Pa, Pf, M, Prf, Prc, Prb, Prab, Prnm, beta, b, f, fab, Tomax, Tmax_ab, MW, eff, y, HVf);

Cp1 = y(4)*(R/MW(4))/(y(4)-1);
Cp2 = y(8)*(R/MW(8))/(y(8)-1);
Tmax = Tomax + CB*(b/bmax)^0.5;
fmax = (1-b)*(1-out2{3}/Tmax)/(eff(4)*HVf/Cp1/Tmax - 1);

g(1) = f - fmax; % f <= fmax
g(2) = (out2{3}+f*HVf/Cp1)/(1+f-b) - Tmax; % Tb <= Tmax
if fab>0
    g(3) = (out2{4} + (f+fab)*HVf/Cp2)/(1+f+fab) - Tmax_ab; % Tab <= Tmaxab
    g(4) = fab - (1+fmax)*(Tmax_ab/out2{4} - 1)/((eff(7)*HVf/Cp2 - Tmax_ab)/out2{4}); % fab <= fmaxab
else
    g(3)=0;
    g(4)=0;
end

ceq = [];
end
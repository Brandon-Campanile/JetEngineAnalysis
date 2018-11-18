function Values = optimization(ST, eType, Nmix, Ta, Pa, Pf, M, ~, ~, Prb, Prab, Prnm, ~, ~, ~, ~, Tomax, Tmax_ab, MW, eff, y, HVf)

A = [];
c = [];
Aeq = [];
beq = [];
T = 0;
TSFC = zeros(1,2);
minVar = cell(1,2);

for ab=0:1
    if eType % tubofan
        if ab % w/ afterburner
            % [bypass ratio, fan pressure ratio, compressor pressure ratio, fuel-air ratio, afterburner fuel-air ratio, bleed ratio]
            lb = [0.001, 1.1, 10, .001, .0005, 0]; % lower bound
            ub = [10, 1.5, 54.545455, .1, .1, .12];        % upper bound
            x0 = [5, 1.25, 40, .05, .03, .05];
        else
            lb = [.001, 1.1, 10, .001, 0, 0];
            ub =[10, 1.5, 54.545455, .1, 0, .12];
            x0 = [5, 1.25, 40, .05, 0, .05];
        end
    else % turbojet
        if ab % w/ afterburner
            lb = [0, 1.1, 10, .001, .0005, 0];
            ub = [0, 1.5, 54.545455, .1, .1, .12];
            x0 = [0, 1.25, 40, .05, .03, .05];
        else
            lb = [0, 1.1, 10, .001, 0, 0];
            ub =[0, 1.5, 54.545455, .1, 0, .12];
            x0 = [0, 1.25, 40, .05, 0, .05];
        end
    end
    
    func = @(x)singleOut(x, T, eType, Nmix, Ta, Pa, Pf, M, Prb, Prab, Prnm, Tomax, Tmax_ab, MW, eff, y, HVf);
    
    nlc = @(x)nonlcon(x, T, ST, eType, Nmix, Ta, Pa, Pf, M, Prb, Prab, Prnm, Tomax, Tmax_ab, MW, eff, y, HVf);
    
    options=optimoptions('fmincon','Display','iter'); % ,'Algorithm','sqp');
    
    b_min = fmincon(func, x0, A, c, Aeq, beq, lb, ub, nlc, options);
    
    minVar{ab+1} = b_min;
    
    TSFC(ab+1) = singleOut(b_min, T, eType, Nmix, Ta, Pa, Pf, M, Prb, Prab, Prnm, Tomax, Tmax_ab, MW, eff, y, HVf);
end

[~,G]=min(TSFC);

Values = minVar{G};

end

function TSFC2 = singleOut(x, T, eType, Nmix, Ta, Pa, Pf, M, Prb, Prab, Prnm, Tomax, Tmax_ab, MW, eff, y, HVf)
beta = x(1);
Prf = x(2);
Prc = x(3);
f = x(4);
fab = x(5);
b = x(6);
output = JetPro_Project(T, eType, Nmix, Ta, Pa, Pf, M, Prf, Prc, Prb, Prab, Prnm, beta, b, f, fab, Tomax, Tmax_ab, MW, eff, y, HVf);
TSFC2 = output(2);
end

function [g,ceq] = nonlcon(x, T, ST, eType, Nmix, Ta, Pa, Pf, M, Prb, Prab, Prnm, Tomax, Tmax_ab, MW, eff, y, HVf)
R=8314;
CB=700;
beta = x(1);
Prf = x(2);
Prc = x(3);
f = x(4);
fab = x(5);
b = x(6);
bmax=.12;

out2 = JetPro_Project(T, eType, Nmix, Ta, Pa, Pf, M, Prf, Prc, Prb, Prab, Prnm, beta, b, f, fab, Tomax, Tmax_ab, MW, eff, y, HVf);

Cp1 = y(4)*(R/MW(4))/(y(4)-1);
Cp2 = y(8)*(R/MW(8))/(y(8)-1);
Tmax = Tomax + CB*(b/bmax)^0.5;
fmax = (1-b)*(1-out2(3)/Tmax)/(eff(4)*HVf/Cp1/Tmax - 1);

g(1) = Prc - 60/Prf; % prc<=60/prf
g(2) = f - fmax; % f <= fmax
g(3) = (out2(3)+f*HVf/Cp1)/(1+f-b) - Tomax - CB*(b/.12)^.5; % Tb <= Tmax
g(4) = ST*1000-out2(1); % ST >= STdesired
if fab>0
    g(5) = (out2(4) + (f+fab)*HVf/Cp2)/(1+f+fab) - Tmax_ab; % Tab <= Tmaxab
    g(6) = fab - (1+fmax)*(Tmax_ab/out2(4) - 1)/((eff(7)*HVf/Cp2 - Tmax_ab)/out2(4)); % fab <= fmaxab
else
    g(5)=0;
    g(6)=0;
end
ceq = [];
end

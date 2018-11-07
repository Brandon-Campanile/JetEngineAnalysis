function Values = optimization(ST, eType, Nmix, Ta, Pa, Pf, M, ~, ~, Prb, Prab, Prnm, ~, ~, ~, ~, Tomax, Tmax_ab, MW, eff, y, HVf)

A = [];
c = [];
Aeq = [];
beq = [];
T = 0;
TSFC = zeros(1,2);
minVar = zeros(1,2);

for ab=0:1
    if eType % tubofan
        if ab % w/ afterburner
            % [bypass ratio, fan pressure ratio, fuel-air ratio, afterburner fuel-air ratio, bleed ratio]
            lb = [0.001, 1.1, .001, .0005, 0]; % lower bound
            ub = [10, 1.5, .1, .1, .5];        % upper bound
            x0 = [5, 1.25, .05, .03, .05];
        else
            lb = [.001, 1.1, .001, 0, 0];
            ub =[10, 1.5, .1, 0, .5];
            x0 = [5, 1.25, .05, 0, .05];
        end
    else % turbojet
        if ab % w/ afterburner
            lb = [0, 1.1, .001, .0005, 0];
            ub = [0, 1.5, .1, .1, .5];
            x0 = [0, 1.25, .05, .03, .05];
        else
            lb = [0, 1.1, .001, 0, 0];
            ub =[0, 1.5, .1, 0, .5];
            x0 = [0, 1.25, .05, 0, .05];
        end
    end
    
    func = @(x)singleOut(x, T, eType, Nmix, Ta, Pa, Pf, M, Prb, Prab, Prnm, Tomax, Tmax_ab, MW, eff, y, HVf);
    
    nlc = @(x)nonlcon(x, T, ST, eType, Nmix, Ta, Pa, Pf, M, Prb, Prab, Prnm, Tomax, Tmax_ab, MW, eff, y, HVf);
    
    b_min = fmincon(func, x0, A, c, Aeq, beq, lb, ub, nlc);
    
    minVar(ab+1) = b_min;
    
    TSFC(ab+1) = singleOut(b_min, T, eType, Nmix, Ta, Pa, Pf, M, Prb, Prab, Prnm, Tomax, Tmax_ab, MW, eff, y, HVf);
end

[~,G]=min(TSFC);

Values = minVar(G);

end

function out = singleOut(x, T, eType, Nmix, Ta, Pa, Pf, M, Prb, Prab, Prnm, Tomax, Tmax_ab, MW, eff, y, HVf)
beta = x(1);
Prf = x(2);
Prc = 60/x(2);
f = x(3);
fab = x(4);
b = x(5);
output = JetPro_Project(T, eType, Nmix, Ta, Pa, Pf, M, Prf, Prc, Prb, Prab, Prnm, beta, b, f, fab, Tomax, Tmax_ab, MW, eff, y, HVf);
out = output(2);
end

function [g,ceq] = nonlcon(x, T, ST, eType, Nmix, Ta, Pa, Pf, M, Prb, Prab, Prnm, Tomax, Tmax_ab, MW, eff, y, HVf)
R=8314;
CB=700;
sbmax=sqrt(.12);
beta = x(1);
Prf = x(2);
Prc = 60/x(2);
f = x(3);
fab = x(4);
b = x(5);
out2 = JetPro_Project(T, eType, Nmix, Ta, Pa, Pf, M, Prf, Prc, Prb, Prab, Prnm, beta, b, f, fab, Tomax, Tmax_ab, MW, eff, y, HVf);
Cp1 = y(4)*(R/MW(4))/(y(4)-1);
g(1) = (HVf*1e6/Cp1 - Tomax)*x(3) + Tomax*x(5) - (CB/sbmax)*sqrt(x(5)) - (CB/sbmax)*x(3)*sqrt(x(5)) + (CB/sbmax)*x(5)*sqrt(x(5)) + (out2(3) - Tomax) ;
g(2) = ST*1000-out2(1);
if fab > 0
    Cp2 = y(8)*(R/MW(8))/(y(8)-1);
    g(3) = out2(4) + (f+fab)*HVf*1e6/Cp2 - (1+f+fab)*Tmax_ab;
end
ceq = [];
end

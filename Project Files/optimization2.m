function output = optimization2(ST, eType, Nmix, Ta, Pa, Pf, M, ~, ~, Prb, Prab, Prnm, ~, ~, ~, ~, Tomax, Tmax_ab, MW, eff, y, HVf)

T = 0;
R = 8314;
lastST = ST;
lastTSFC = 1;

for beta = 7:.5:10
    fprintf('beta = %d\n',beta);
    for Prf = 1.35:.025:1.5
        Prc = 60/Prf;
        for f=.0005:.005:.035
            for fab=0.0005:.005:.035
                for b=0.08:.006:.11
                    fprintf('b = %d\n',b/.03);
                    out = JetPro_Project(T, eType, Nmix, Ta, Pa, Pf, M, Prf, Prc, Prb, Prab, Prnm, beta, b, f, fab, Tomax, Tmax_ab, MW, eff, y, HVf);
                    CB = 700;
                    
                    Cp1 = y(4)*(R/MW(4))/(y(4)-1);
                    Cp2 = y(8)*(R/MW(8))/(y(8)-1);
                    
                    Temp = (out(3)+f*HVf/Cp1)/(1+f-b);
                    Tmax = Tomax + CB*(b/.12)^.5;
                    Tab = (out(4)+(f+fab)*HVf/Cp2)/(1+f+fab);
                    
                    fmax = (1-b)*(1-out(3)/Tomax)/(eff(4)*HVf/Cp1/Tomax - 1);
                    fmax_ab = (1+fmax)*(Tmax_ab/out(4) - 1)/((eff(7)*HVf/Cp2 - Tmax_ab)/out(4));
                    
                    if Temp<=Tmax && f<=fmax && Tab<=Tmax_ab && fab<=fmax_ab && out(1)>=lastST && 0<out(2) && out(2)<lastTSFC
                        values = [beta, Prf, f, fab, b, out(1), out(2)];
                        lastST = out(1);
                        lastTSFC = out(2);
                    end
                end
            end
        end
    end
end
if exist('values', 'var')
    output = values;
else
    error('could not converge')
end
end

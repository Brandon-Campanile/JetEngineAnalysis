function b_min = optimization(ST, eType, Nmix, Ta, Pa, Pf, M, ~, ~, Prb, Prab, Prnm, ~, ~, ~, ~, Tomax, Tmax_ab, MW, eff, y, HVf)

    % [bypass ratio, fan pressure ratio, fuel-air ratio, afterburner fuel-air ratio, bleed ratio]
    lb = [0 1.1 0 0 0];         % lower bound 
    ub = [15 1.5 .2 .2 1];      % upper bound
    A = [];
    c = [];
    Aeq = [];
    beq = [];
    x0 = [10 1.3 .04 .015 .02];
    T = 0;
    
    func = @(x)singleOut(x, T, eType, Nmix, Ta, Pa, Pf, M, Prb, Prab, Prnm, Tomax, Tmax_ab, MW, eff, y, HVf);
    
    nlc = @(x)nonlcon(x, T, ST, eType, Nmix, Ta, Pa, Pf, M, Prb, Prab, Prnm, Tomax, Tmax_ab, MW, eff, y, HVf);
    
%    options = optimoptions(@fmincon, 'OutputFcn', @outfun);
 
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
    
%     function stop = outfun(func, ~, ~)
%         beta = func(1);
%         Prf = func(2);
%         Prc = 60/func(2);
%         f = func(3);
%         fab = func(4);
%         b = func(5);
%         
%         this = JetPro_Project(T, eType, Nmix, Ta, Pa, Pf, M, Prf, Prc, Prb, Prab, Prnm, beta, b, f, fab, Tomax, Tmax_ab, MW, eff, y, HVf);
%         iST = this(1);
%         if iST<=ST
%             stop= true;
%         else
%             stop= false;
%         end
%     end

    

    function [c,ceq] = nonlcon(x, T, ST, eType, Nmix, Ta, Pa, Pf, M, Prb, Prab, Prnm, Tomax, Tmax_ab, MW, eff, y, HVf)
        R=8314;
        CB=700;
        sbmax=sqrt(.12);
        beta = x(1);
        Prf = x(2);
        Prc = 60/x(2);
        f = x(3);
        fab = x(4);
        b = x(5);
        
        out = JetPro_Project(T, eType, Nmix, Ta, Pa, Pf, M, Prf, Prc, Prb, Prab, Prnm, beta, b, f, fab, Tomax, Tmax_ab, MW, eff, y, HVf);
        Cp1 = y(4)*(R/MW(4))/(y(4)-1);
        c(1) = (HVf/Cp1 - Tomax)*x(3) + Tomax*x(5) - (CB/sbmax)*sqrt(x(5)) - (CB/sbmax)*x(3)*sqrt(x(5)) + (CB/sbmax)*x(5)*sqrt(x(5)) + (out(3) - Tomax) ; 
        
        c(2) = ST*1000-out(1);
        
        if fab > 0
            Cp2 = y(8)*(R/MW(8))/(y(8)-1);
            c(3) = out(4) + (f+fab)*HVf/Cp2 - (1+f+fab)*Tmax_ab;
            ceq = [];
        end
    end

    b_min = fmincon(func, x0, A, c, Aeq, beq, lb, ub, nlc);

end
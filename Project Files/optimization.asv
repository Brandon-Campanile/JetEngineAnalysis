function b_min = optimization(ST, eType, Nmix, Ta, Pa, Pf, M, ~, ~, Prb, Prab, Prnm, ~, ~, ~, ~, Tomax, Tmax_ab, MW, eff, y, HVf)

    % [bypass ratio, fan pressure ratio, fuel-air ratio, afterburner fuel-air ratio, bleed ratio]
    lb=[0 1.1 0 0 0];       % lower bound 
    ub=[15 1.5 .2 .2 1];      % upper bound
    A=[];
    c=[];
    Aeq=[];
    beq=[];
    x0=[10 1.3 .04 .03 .05];
    
    function out = singleOut(eType, Nmix, Ta, Pa, Pf, M, Prf, Prc, Prb, Prab, Prnm, beta, b, f, fab, Tomax, Tmax_ab, MW, eff, y, HVf)
        output = JetPro_Project(eType, Nmix, Ta, Pa, Pf, M, Prf, Prc, Prb, Prab, Prnm, beta, b, f, fab, Tomax, Tmax_ab, MW, eff, y, HVf);
        out = output(2);
    end
    
    func = @(x) singleOut(eType, Nmix, Ta, Pa, Pf, M, x(2), 60/x(2), Prb, Prab, Prnm, x(1), x(5), x(3), x(4), Tomax, Tmax_ab, MW, eff, y, HVf);
    
    function stop = outfun(func, ~, ~)
        this = JetPro_Project(eType, Nmix, Ta, Pa, Pf, M, func(2), 60/func(2), Prb, Prab, Prnm, func(1), func(5), func(3), func(4), Tomax, Tmax_ab, MW, eff, y, HVf);
        iST = this(1);
        if iST<=ST
            stop= true;
        else
            stop= false;
        end
    end

    function [c,ceq] = nonlcon1(x)
        R=8314;
        CB=700;
        sbmax=sqrt(.12);
        out = JetPro_Project(eType, Nmix, Ta, Pa, Pf, M, x(2), 60/x(2), Prb, Prab, Prnm, x(1), x(5), x(3), x(4), Tomax, Tmax_ab, MW, eff, y, HVf);
        Cp = y(4)*(R/MW(4))/(y(4)-1);
        c = (HVf/Cp - Tomax)*x(3) + Tomax*x(5) - (CB/sbmax)*sqrt(x(5)) - (CB/sbmax)*x(3)*sqrt(x(5)) + (CB/sbmax)*x(5)*sqrt(x(5)) + (out(3) - Tomax) ; 
        ceq = [];
    end

    function [c2,ceq2] = nonlcon2(x2)
        R=8314;
        CB=700;
        sbmax=sqrt(.12);
        out = JetPro_Project(eType, Nmix, Ta, Pa, Pf, M, x2(2), 60/x2(2), Prb, Prab, Prnm, x2(1), x2(5), x2(3), x2(4), Tomax, Tmax_ab, MW, eff, y, HVf);
        Cp = y(4)*(R/MW(4))/(y(4)-1);
        c2 = (HVf/Cp - Tomax)*x2(3) + Tomax*x2(5) - (CB/sbmax)*sqrt(x2(5)) - (CB/sbmax)*x2(3)*sqrt(x2(5)) + (CB/sbmax)*x2(5)*sqrt(x2(5)) + (out(3) - Tomax) ; 
        ceq2 = [];
    end

    options = optimoptions(@fmincon, 'OutputFcn', @outfun);

    b_min = fmincon(func, x0, A, c, Aeq, beq, lb, ub, {@nonlcon1, @nonlcon2}, options);

end
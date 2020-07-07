function [p, fminres, fitresult] = My_two_gaussianfit_constant_ratio(xdata, ydata, varargin)

    ratio=1.6757;
    % obtain min and max value for determining fitting initial para.
    miny = min(ydata(:));
    maxy = max(ydata(:));
    
    opts = optimset('Display','off');
        
    % fitting for x axis projection
    if ~isempty(length(varargin)) 
        init_guess = varargin{1}; 
    else
        init_guess = [miny maxy-miny mean(xdata)*1/4 length(xdata)/8 miny maxy-miny mean(xdata)*3/4 length(xdata)/8];
    end
    init_guess=init_guess([1:4 6]);
    lb=init_guess*0.1;
    ub=init_guess*10;
    fun = @(p,xdata) abs(p(1))*exp(-((xdata-p(2))/p(3)).^2) + abs(p(4))*exp(-((xdata-p(2)*ratio)/p(5)).^2);
    [p,fminres] = lsqcurvefit(fun,init_guess,xdata,ydata,lb,ub,opts);
    p(6)=p(5);
    p(5)=p(2)*ratio;
    fitresult = abs(p(1))*exp(-((xdata-min(p(2),p(5)))/p(3)).^2) + abs(p(4))*exp(-((xdata-max(p(2),p(5)))/p(6)).^2);
    
end
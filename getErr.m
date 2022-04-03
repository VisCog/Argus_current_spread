function [err,x,prob] =  getErr(p,dte,a,dta,resp,funName)
% [err,x,prob] =  getErr(p,a,d,resp,funName)
% funName: 'normcdf', 'logit', or 'weibull'

if ~exist('funName')
    funName = 'logit';
end

switch funName
    case 'normcdf'
        if ~isfield(p,'sd')
            p.sd = 1;
        end
        x= p.b0 + p.ka*a+p.kdta*dta+p.kdte*dte;
        prob = normcdf(x,0,p.sd);
    case 'logit'
        x= p.b0 + p.ka*a+p.kdta*dta+p.kdte*dte;
        prob = 1./(1+exp(-x));
end

% prob = max(min(prob, 0.995), 0.005);
if ~isempty(resp) 
    % calculate negative maximum log likelihood
    err = -sum(resp.*log(prob) + (1-resp).*log(1-prob));
else
    err = NaN;
end


function [rho, pval, uniform] = circ_corrcc_uniform(alpha1, alpha2, test_uniform, test_uniform_alpha)
%
% [rho, pHodgesAjne] = circ_corrcc_uniform(alpha1, alpha2, test_uniform, test_uniform_alpha)
%   Circular correlation coefficient for two periodic signals with
%   uniformly distributed phase directions
%
%   Input:
%     alpha1	            sample of angles in radians
%     alpha2                sample of angles in radians
%     test_uniform          perform Hodges-Ajne test for uniform distribution of input signals (default = 1)
%     test_uniform_alpha    alpha threshold for Hodges-Ajne test (default = .05)
%
%   Output:
%     rho                   correlation coefficient
%     pval                  p value for correlation
%     uniform               structure with Hodges-Ajne test (circ_otest) for non-uniformity of input signals
%        .h                 test result
%        .p_alpha1          p value alpha 1
%        .p_alpha2          p value alpha 2
%
% References:
%   Circular correlation in EEG hyperscanning:
%       Zimmermann et al. (under revision). Arbitrary methodological decisions skew inter-brain synchronization estimates in hyperscanning-EEG studies.
%   CircStat toolbox:
%       Berens, P. (2009). CircStat: A MATLAB Toolbox for Circular Statistics. Journal of Statistical Software, 31(10), 1–21. https://doi.org/10.18637/jss.v031.i10
%   Topics on circular statistics:
%       Topics in circular statistics, S.R. Jammalamadaka et al., Chapter 8.ii (formula 8.2.4) 
%

% adapted from Circular Statistics Toolbox for Matlab
%   by Philipp Berens, 2009
%   berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
%
% edited by Marius Zimmermann, 2024, University of Regensburg
%   marius.zimmermann@ur.de
%
% Versions
%   PHB   6/7/2008
%   MZ    20/9/2024 (adjusted for uniform distributions, e.g. EEG data)
%

%% check input arguments and default values
arguments
    alpha1  double {mustBeVector}
    alpha2  double {mustBeVector,mustBeEqualSize(alpha1,alpha2)}
    test_uniform {mustBeNumeric,mustBeNonempty,mustBeMember(test_uniform,[0 1])} = 1;
    test_uniform_alpha {mustBeNumeric,mustBeNonempty,mustBeInRange(test_uniform_alpha,0,1,'exclude-lower')} = .05;
end

if ~iscolumn(alpha1)
    alpha1 = reshape(alpha1,[length(alpha1) 1]);
    alpha2 = reshape(alpha2,[length(alpha2) 1]);
end

%% perform Hedges-Ajne test for non-uniformity of input data (using circ_otest)
% H0: the population is uniformly distributed around the circle
% HA: the population is not uniformly distributed around the circle

if test_uniform
    pHodgesAjne = nan(1,2);
    pHodgesAjne(1) = circ_otest(alpha1);
    pHodgesAjne(2) = circ_otest(alpha2);

    uniform.h = any(pHodgesAjne < test_uniform_alpha);
    uniform.p_alpha1 = pHodgesAjne(1);
    uniform.p_alpha2 = pHodgesAjne(2);
    
    if any(pHodgesAjne < test_uniform_alpha)
        warning(sprintf('input signal(s) not uniformly distributed! p[alpha1] = %.3f, p[alpha2] = %.3f',pHodgesAjne))
    end
else
    uniform = NaN;
end


%% calc circular correlations for uniform distrubution
% using:
% m - n = circ_mean(wrapToPi(alpha1-alpha2));
% m + n = circ_mean(wrapToPi(alpha1+alpha2));
% thus ...
% -2n = circ_mean(wrapToPi(alpha1-alpha2)) - circ_mean(wrapToPi(alpha1+alpha2));

n = -1*(circ_mean(wrapToPi(alpha1-alpha2)) - circ_mean(wrapToPi(alpha1+alpha2)))/2;
m = circ_mean(wrapToPi(alpha1-alpha2))+n;

x_sin = sin(alpha1 - m);
y_sin = sin(alpha2 - n);

r_minus = abs(sum(exp((alpha1-alpha2) * 1i)));
r_plus  = abs(sum(exp((alpha1+alpha2) * 1i)));

num     = (r_minus - r_plus);
den     = 2*sqrt( sum(x_sin .^2) .* sum(y_sin .^2));

rho = num / den;

%% compute p value
l20 = mean(sin(alpha1 - circ_mean(wrapToPi(alpha1))).^2);
l02 = mean(sin(alpha2 - circ_mean(wrapToPi(alpha2))).^2);
l22 = mean((sin(alpha1 - circ_mean(wrapToPi(alpha1))).^2) .* (sin(alpha2 - circ_mean(wrapToPi(alpha2))).^2));

ts = sqrt((length(alpha1) * l20 * l02)/l22) * rho;
pval = 2 * (1 - normcdf(abs(ts)));

end


%% Custom validation function (source: matlab help center)
function mustBeEqualSize(a,b)
    % Test for equal size
    if ~isequal(size(a),size(b))
        eid = 'Size:notEqual';
        msg = 'Size of first input must equal size of second input.';
        error(eid,msg)
    end
end

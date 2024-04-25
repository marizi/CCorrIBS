function [rho] = circ_corrcc_uniform(alpha1, alpha2)
%
% [rho pval ts] = circ_corrcc_uniform(alpha1, alpha2)
%   Circular correlation coefficient for two periodic signals with
%   uniformly distributed phase directions
%
%   Input:
%     alpha1	sample of angles in radians
%     alpha2	sample of angles in radians
%
%   Output:
%     rho     correlation coefficient
%     pval    p-value
%
% References:
%   Topics in circular statistics, S.R. Jammalamadaka et al., Chapter 8.ii (formula 8.2.4) 
%
% PHB   6/7/2008
% MZ    1/8/2020 (uniform)
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html
%
% edited by Marius Zimmermann, 2020 (marz@dtu.dk)


if size(alpha1,2) > size(alpha1,1)
	alpha1 = alpha1';
end

if size(alpha2,2) > size(alpha2,1)
	alpha2 = alpha2';
end

if length(alpha1)~=length(alpha2)
  error('Input dimensions do not match.')
end

% calc circular correlations for uniform distrubution

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

% % compute pvalue
% n = length(alpha1);
% alpha1_bar = circ_mean(alpha1);
% alpha2_bar = circ_mean(alpha2);
% l20 = mean(sin(alpha1 - alpha1_bar).^2);
% l02 = mean(sin(alpha2 - alpha2_bar).^2);
% l22 = mean((sin(alpha1 - alpha1_bar).^2) .* (sin(alpha2 - alpha2_bar).^2));

% ts = sqrt((n * l20 * l02)/l22) * rho;
% pval = 2 * (1 - normcdf(abs(ts)));


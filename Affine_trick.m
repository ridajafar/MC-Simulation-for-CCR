function [A_vec, B_vec] = Affine_trick(t, pricing_grid, a, sigma, ZC_curve )
% Computes A and B vectors
%
% INPUT
% t:            simulation grid for the MC simulation
% pricing_grid: time grid of the pricing dates
% a:            Rates Mean reversion speed parameter
% sigma:        Short-rate volatility 'sigma' on the slides
% ZC_curve:     bootstrapped ZC curve


% Set time parameters
index = find([0 pricing_grid]==t);
dT = pricing_grid(index:end) - t;
pricing_grid = [0 pricing_grid];
t_shift = t + 0.0001;

% Compute shifted and non shifted ZC curves
ZC_pricing = [ZC_curve(1,2); ZC_curve(1:length(pricing_grid)-1,2)];
ZC_shift = interp1(pricing_grid,ZC_pricing, t_shift, 'linear');

% Compute forward
f = (ZC_shift*t_shift - ZC_pricing(index)*t)./(t_shift-t);

% Compute A_vec and B_vec
A_vec = zeros(1,length(pricing_grid)-1);
B_vec = zeros(1,length(pricing_grid)-1);
B_vec(index:end) = (1-exp(-a*dT))/a;
A_vec(index:end) = exp(-ZC_pricing(index+1:end)'.*pricing_grid(index+1:end) + ZC_pricing(index)*t) .* exp(B_vec(index:end)*f - sigma^2/(4*a) * (1-exp(-2*a*t)) * B_vec(index:end).^2) ;

end
function [ mean_add, mean_mult, sq_var, dt ] = mc_param( sim_grid, ZC_curve, a, sigma )
% Computes MC parameters
%
% INPUT
% sim_grid: simulation grid for the MC simulation
% ZC_curve: bootstrapped ZC curve
% a:        Rates Mean reversion speed parameter
% sigma:    Short-rate volatility 'sigma' on the slides


% Set shifted and non shifted time vectors
sim_grid = [0 sim_grid];
sim_grid_shift = sim_grid + 0.0001;
dt = sim_grid(2:end) - sim_grid(1:end-1);

% Compute shifted and non shifted ZC curves
ZC_sim = [ZC_curve(1,2), ZC_curve(1:length(sim_grid)-1,2)'];
ZC_shift = interp1(sim_grid, ZC_sim, sim_grid_shift, 'linear', 'extrap');

% Compute shifted and non shifted Discounts
B_shift = exp(-ZC_shift.*sim_grid_shift);
B = exp(-ZC_sim.*sim_grid);

% Compute Forwards
f = -(log(B_shift./B))./(sim_grid_shift - sim_grid);

% Compute alpha
alpha = f + sigma^2/(2*a^2)*(1-exp(-a*sim_grid)).^2;

% Compute outputs
mean_add = alpha(2:end) - alpha(1:end-1).*exp(-a*dt);
mean_mult = exp(-a*dt);
sq_var = sqrt(sigma^2/2/a * (1 - exp(-2*a*dt)));

end
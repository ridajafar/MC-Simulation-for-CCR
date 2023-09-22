%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Template Case RM 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Montecarlo simulation of the interest rate term structure 
%   for CCR purposes
%   Rate dynamics follows the Hull-White model
%   (See Brigo-Mercurio 3.33 page 73)
%
%   INPUT Market data are stored as:
%
%   1. Market_data: table of ICAP quotes for the strip of IRS
%       Column #1: IRS maturity (year frac)
%       Column #2: MID rate
%
%   2. sim_grid: Grid for simulation (year fractions) - t_i on the slides
%   
%   3. a:      Rates Mean reversion speed parameter 'a' on the slides
%   4. sigma:  Short-rate volatility 'sigma' on the slides
%
%   5. IRS_cf:  Table of cashflows of the "par yield bond" used to simulate
%       Column #1: cashflow maturity (year frac)
%       Column #2: cashflow amount
%
%    OUTPUT data are stored as:
%   6. ZC_curve: Table of ZC rates (cont. comp. 30/360)
%      Maturities are year fractions
%   7. A set of scalar/vector variables with self-explanatory names
%   
%   8. Required functions' template:
%   8.1. ZC_bootstrap_IRS_only: Analytical bootstrap of ZC curve from IRS
%   function [ZC_curve]=ZC_bootstrap_IRS_only(market_data, frequency)
%
%   8.2. mc_param: (time saving) pre-calculated parameters needed to build the conditional
%   short-term rate r(t) given a realization r(s) (s < t) on the sim_grid.
%   Given r(s) and a random realization z under N(0,1)then 
%   r(t) = r(s) * mean_mult + mean_add + sqrt(variance) * z
%   See Brigo-Mercurio 3.37  
%   The four output of this function are vectors with dimension equal to
%   the number of points of the time grid (i.e. same dimension as the input
%   vector sim_grid). 
%   On an equally-spaced time grid some parameters are constant across
%   time.
%   function [ mean_add, mean_mult, sq_var, dt ] = mc_param( sim_grid, ZC_curve, a, sigma )
%
%   8.3 Affine_trick: 
%   Pricing of a ZCB issued at t with maturity T is given by Brigo-Mercurio
%   3.39: P(t,T) = A(t,T) * exp(-B(t,T)*r(t)) 
%   Given t and a vector of future maturities T_j (stored in input vector
%   pricing_grid), corresponding values of A(t,T) and B(t,T) are stored in 
%   the (vector) function outputs: A_vec, B_vec
%   function [ A_vec, B_vec ] = Affine_trick( t, pricing_grid, a, sigma, ZC_curve )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

%   Rates Mean reversion speed parameter 'a' on the slides
a = 0.0195;   

%   Short-rate volatility 'sigma' on the slides
sigma = 0.0086;
% sigma = 0.0172; %Classroom discussion

% IRS Quotes
market_data = [ 0.25 4.84; 0.50 5.10; 0.75 5.27; 1.00 5.43; 
                1.25 5.32; 1.50 5.22; 1.75 5.11; 2.00 5.00; 
                2.25 4.90; 2.50 4.80; 2.75 4.70; 3.00 4.60; 
                3.25 4.54; 3.50 4.48; 3.75 4.41; 4.00 4.35; 
                4.25 4.31; 4.50 4.27; 4.75 4.23; 5.00 4.20; 
                5.25 4.17; 5.50 4.15; 5.75 4.13; 6.00 4.11; 
                6.25 4.09; 6.50 4.07; 6.75 4.05; 7.00 4.03; 
                7.25 4.02; 7.50 4.01; 7.75 4.00; 8.00 4.00; 
                8.25 3.97; 8.50 3.97; 8.75 3.96; 9.00 3.95; 
                9.25 3.94; 9.50 3.94; 9.75 3.93; 10.0 3.92];


% IRS cashflows (fixed leg only, and treated as a fixed rate bond)
T_max = 8.0;  % IRS Maturity
j_max = 32;    % Number of future cashflows of the IRS (only considering fixed leg)
cpn = 0.01;  % Coupon amount
%cpn = 0.02;   % Coupon amount Q7
irs_fixed_coupon_freq = 4;  %Annual fixed coupons 
% j_max = T_max * irs_fixed_coupon_freq;    % Number of future cashflows of the IRS (only considering fixed leg)
% cpn = 0.01;  % Coupon amount
IRS_cf(:,1) = linspace(0.25,T_max,j_max);
IRS_cf(:,2) = -cpn;
IRS_cf(j_max,2) = IRS_cf(j_max,2) - 1.0; % fixed leg treated as a fixed rate bond
% (workplace) vector (row) of time points of each IRS cashflow
pricing_grid = IRS_cf(:,1)'; % Let us use te symbol T_j for each point of the pricing grid

% Vector (row) tith time-grid for simulation (year fractions) - t_i on the slides
i_max = j_max; % Number of grid points
sim_grid = linspace(0.25,T_max,i_max);
EE = zeros(1,i_max);      % Expected Exposure on the simulation gridpoints
PFE = zeros(1,i_max);     % Pot. Fut. Exposure on the simulation gridpoints

% Number of mc iterations
N = 250000;
rng(42)

% Confidence level of PFE
alpha = 0.95;
k = round((1-alpha)*N); % index of the sorted MtF corresponding to alpha

%% Bootstrap
[ZC_curve]=ZC_bootstrap_IRS_only(market_data,irs_fixed_coupon_freq);

%%
% pre-calculated parameters
[mean_add, mean_mult, sq_var, dt] = mc_param( sim_grid, ZC_curve, a, sigma );

%%
% Q1 Test on the initial MtM of the IRS
% Also a way to make some exercise on the usage of the "Affine trick"
% P_vec: Pricing of a ZCB issued at t with maturity T_j i.e. P(t,T_j)
t = 0.0;
spot_rate = ZC_curve(1,2);
[A_vec, B_vec] = Affine_trick(t, pricing_grid, a, sigma, ZC_curve );
P_vec = A_vec.*exp(-B_vec*spot_rate);   % Vector of discount factors, one element for each future cashflow
t0_MtM = sum(P_vec' .* IRS_cf(:,2))+1; % Payer
%t0_MtM = -sum(P_vec' .* IRS_cf(:,2))-1; % Receiver
disp('––– Q1: Current MtM of the IRS –––')
fprintf('Current MtM of the IRS: %.2f \n', t0_MtM)
disp(' ')

%%
% Montecarlo: N scenarios, composed bt j_max time steps each
MtF = zeros(1,N);       % Mark-to-future
r = zeros(1,N);         % Simulated short-term rate

% Iteration on the steps of the simulation grid
for i = 1:i_max
    
    % Iteration on the steps of the simulation grid
    t = sim_grid(i);
    
    % Generate N standard 1-d normal variables stored to z(1,N)
    z = randn(1,N);
    
  % Scenarios for the short-term rate on the simulation grid 
  % stored to r(1,N)
    if i ==1
        r = spot_rate*mean_mult(i) + mean_add(i) + sq_var(i)*z; %Insert here the rate simulated for the first time step. All scenarios take spot_rate as baseline
    else
        r = r*mean_mult(i) + mean_add(i) + sq_var(i)*z ; %Insert here the rate simulated for the other time styeps. Each scenario takes r at the previous time step as baseline
    end
    
    % Given the N simulated short-term rates on the simulation grid for
    % the j-th path of the Montecarlo simulation, stored in r,
    % the affine structure of the HW model allows the analytical calculation of the simulated MtM 
    % of the portfolio at the i-th step of the simul. grid (t_i stored in
    % the variable t)
    [A_vec, B_vec] = Affine_trick(t, pricing_grid, a, sigma, ZC_curve );
    % Be careful! As t_i increases, if T_j <= t_i, the corresponding 
    % element of A_vec and of B_vec must be set to 0.
    % Indeed as the time step of the simulation increaes, there are 
    % fewer cashflows which matter for the MtF
    
    % Vector of simulated MtM of the portfolio at the i-th step of the simul. grid
    P_vec = @(r) A_vec.*exp(-B_vec*r);
    MtF = arrayfun(@(r) sum(P_vec(r)' .* IRS_cf(:,2))+1*(t<T_max),r); %Payer Fix
    %MtF = arrayfun(@(r) -sum(P_vec(r)' .* IRS_cf(:,2))-1*(t<T_max),r); %Receiver Fix

    % Risk measures:

    % Q2. Expected Exposure
    if max(MtF)>0
        EE(i) = mean(MtF(MtF>0));
    end
    
    % Q3. Potential Future Exposure
    MtF_ordered = sort(MtF,'descend');
    PFE(i)=MtF_ordered(k); % 

    % Once the risk measures have been derived at t_i, we trash MtF to spare
    % memory
end

%% Summary

% Q4 Expected Positive Exposure
EPE=mean(EE);
disp('––– Q4: EPE of the IRS –––')
fprintf('EPE as a percentage of the notional: %.1f%%\n', 100*EPE)
disp(' ')

% Q5 Peak - PFE
Peak_PFE = max(PFE);
disp('––– Q5: Peak PFE of the IRS –––')
fprintf('Peak PF as a percentage of the notional: %.1f%%\n', 100*Peak_PFE)
disp(' ')


%% Chart
figure
plot(EE,'o','LineWidth',2)
hold on
plot(PFE,'*','LineWidth',2,'MarkerEdgeColor','r')
hold on
plot(sim_grid*4,EPE,'+','LineWidth',2,'MarkerEdgeColor','g')
xlim([0 length(EE)+1])
% Set limits to the vertical axis. Should the two charts be identical?
% ylim([0 xx]) % Receiver fix
%ylim([0 yy]) % Payer fix
legend('Expected Exposure','Potential Future Exposure','Expected Positive Exposure')
%legend('Potential Future Exposure (PFE)','Peak PFE')
xlabel('time (quarters)')
title('CCR Risk Measures for 10y IRS - Co-terminal calibration')
set(gca,'ygrid','on')
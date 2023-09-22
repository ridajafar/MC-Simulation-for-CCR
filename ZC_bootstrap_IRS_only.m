function [ZC_curve]=ZC_bootstrap_IRS_only(IRS_data, freq)
% Bootstrap of ZC curve from IRS
% 
% INPUT:
% Market_data:  table of ICAP quotes for the strip of IRS
%   Column #1:  IRS maturity (year frac)
%   Column #2:  MID rate
% freq:         frequency of coupons
% 
% OUTPUT:
% ZC_curve:     Table of ZC rates (cont. comp. 30/360)
%               Maturities are year fractions

% Swap rates
swap_rates = IRS_data(:,2)/100;

% Compute discount factors from swap rates
B = zeros(size(IRS_data,1),1);
for i=1:size(IRS_data,1)
    B(i) = (1 - swap_rates(i)*sum(B(1:i-1))/freq)/(1 + swap_rates(i)/freq);
end

% Compute zero rates from discount factors
ZC_curve(:,1) = IRS_data(:,1); 
ZC_curve(:,2) = (-log(B))./ZC_curve(:,1);

end
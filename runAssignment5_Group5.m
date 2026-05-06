% runAssignment5
% group 5, AY2025-2026
%
% Calibration of spot volatilities on LMM, delta, vega sensitivities of a
% structured bond, delta and vega hedging by using swaps and caps, pricing
% an exoticc cup on a BMM

addpath('bootstrap',genpath('ex_1'));
formatData='dd/mm/yyyy';
[datesSet, ratesSet] = readExcelData_MAC('MktData_CurveBootstrap.xls', formatData);
[dates, discounts, ~] = bootstrap(datesSet, ratesSet); 
% Load the Eurostoxx market data from the .mat file
data = load('eurostoxx_Poli.mat'); 
data = data.cSelect;
strike = 3200;
spread = 0.013;
notional = 1e8;
parameters = compute_parameters(data,dates,discounts,strike,spread);


%% Case study 1

sigma = 0.2;
k = 1;
eta = 3;
N_sim = 1e7;
alpha = 0.5;
first_coupon = 0.06;
second_coupon = 0.02;

upfront = compute_upfront(parameters,sigma,k,eta,...
    N_sim,alpha,first_coupon,second_coupon,spread,notional);
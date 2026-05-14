% runAssignment6
% group 5, AY2025-2026
%
% Pricing of a 2y/3y EURO STOXX certificate swap under a calibrated NIG model,
% with Black smile-adjusted and flat benchmarks; pricing of a 10y Bermudan
% Payer Swaption (non-call 2) under a 1-factor Hull-White trinomial tree.

addpath('bootstrap',genpath('ex_1'),genpath('ex_2'));
formatData='dd/mm/yyyy';
[datesSet, ratesSet] = readExcelData_MAC('MktData_CurveBootstrap.xls', formatData);
[dates, discounts, ~] = bootstrap(datesSet, ratesSet); 
% Load the Eurostoxx market data from the .mat file
data = load('eurostoxx_Poli.mat'); 
data = data.cSelect;
strike = 3200;
spread = 0.013;
notional = 1e8;
volatilities = data.surface(:);
strikes = data.strikes(:);
[parameters,sigma,k,eta] = compute_parameters(data,dates,discounts, ...
    strike,spread,volatilities,strikes);

%% Case study 1

N_sim = 1e7;
alpha = 0.5;
first_coupon = 0.06;
second_coupon = 0.02;

% Point a)
upfront_nig = compute_upfront(parameters,sigma,k,eta,...
    N_sim,alpha,first_coupon,second_coupon,spread,notional, 'NIG');

% Point b)
upfront_black_smile = compute_upfront(parameters,sigma,k,eta,...
    N_sim,alpha,first_coupon,second_coupon,spread,notional,'BlackSmile');

% Point d)
upfront_3y_nig = compute_upfront_3y(parameters, sigma, k, eta, ...
    N_sim, alpha, first_coupon, second_coupon, spread, notional);

% Point e)
upfront_black_flat_2y = compute_upfront(parameters, sigma, k, eta, ...
    N_sim, alpha, first_coupon, second_coupon, spread, notional, 'BlackFlat');

fprintf('\n--- Case Study 1: 2Y Certificate Swap ---\n');
fprintf('Upfront 2Y_NIG         = %.8f\n', upfront_nig);
fprintf('Upfront Black Smile = %.8f\n', upfront_black_smile);
fprintf('Upfront 3Y_NIG         = %.8f\n', upfront_3y_nig);
fprintf('Upfront Black Flat  = %.8f\n', upfront_black_flat_2y);

fprintf('Black Flat  - NIG = %.8f  (%.4f bps)\n', ...
    upfront_black_flat_2y - upfront_nig, ...
    (upfront_black_flat_2y - upfront_nig) * 1e4);

%% Case Study 2: Bermudan Swaption Pricing via Hull-White Model

settlement_date = datetime(2008, 2, 15);
% Calculate the start date as Settlement + 2 business days
startDate_num = business_date_offset(settlement_date, 'day_offset', 2);
startDate = datetime(startDate_num, 'ConvertFrom', 'datenum');

% Hull-White Model parameters
a = 0.11;       % Mean reversion speed
sigma = 0.008;  % Volatility of the short rate
K = 0.05;       % Fixed strike rate (coupon) of the underlying swap

% POINT A: PRICING

% Convergence settings: test different number of time steps per year
precision_levels = [1, 4, 12, 52, 365]; 

% Generate Exercise dates: annually from Year 2 to Year 9
y_ex = (2:9)';
exercise_dates = datetime(arrayfun(@(y) business_date_offset(startDate, ...
    'year_offset', y), y_ex), 'ConvertFrom', 'datenum');

% Generate Underlying Swap payment dates: annually from Year 1 to Year 10
y_pay = (1:10)';
swap_payment_dates = datetime(arrayfun(@(y) business_date_offset(startDate,...
    'year_offset', y), y_pay), 'ConvertFrom', 'datenum');

% We compute the prices with different precision level
[bermudan_prices, summary] = run_hw_pricing_swaption(a, sigma, K, startDate,...
    exercise_dates, swap_payment_dates, precision_levels, dates, discounts, 'Bermudan');

% Display the results
disp(summary);

% POINT B: CHECK CORRECTED TREE IMPLEMENTATION

disp(' POINT B: TREE CALIBRATION & SANITY CHECKS');

disp('1. Test 10Y ZCB Calibration (Tree Convergence vs Market Curve)');

% 10-Year Maturity (which is the last date in swap_payment_dates)
maturity_zcb = swap_payment_dates(end); 

% Call the master function with the 'ZCB' flag and empty exercise dates
[zcb_tree_prices, zcb_table] = run_hw_pricing_swaption(a, sigma, K, startDate,...
    datetime([], 'ConvertFrom', 'datenum'), maturity_zcb, precision_levels, dates, discounts, 'ZCB');

% Compute the exact analytical target from the initial market curve
market_price_zcb = get_discount_factor_by_zero_rates_linear_interp(...
    datenum(startDate), datenum(maturity_zcb), datenum(dates), discounts);

% Display the ZCB convergence table and final error
disp(zcb_table);
fprintf('Target Market Curve Price: %.6f\n', market_price_zcb);
fprintf('Error at 365 steps: %.2e\n\n', abs(zcb_tree_prices(end) - market_price_zcb));


% CHECK 2: EUROPEAN SWAPTION (9Y) 
% We test a long-dated European Swaption. This forces the tree to perform
% stochastic backward induction over a 9-year horizon, deeply testing 
% the tree's geometry and probabilities.
disp('2. Test 9Y European Swaption (Tree Convergence vs Jamshidian)');

% Isolate the last exercise date (Year 9) to force European style
eur_ex_date = exercise_dates(end); 

% Filter the swap payment dates to include only those after the exercise date
% For a 9Y exercise, this will only be the final payment at Year 10
eur_swap_payments = swap_payment_dates(swap_payment_dates > eur_ex_date);

% Call the master function with the 'European' flag and a single exercise date
[eur_tree_prices, eur_table] = run_hw_pricing_swaption(a, sigma, K, startDate,...
    eur_ex_date, eur_swap_payments, precision_levels, dates, discounts, 'European');

% Compute the exact analytical target using Jamshidian's decomposition
eur_analytical_price = price_jamshidian(a, sigma, K, eur_ex_date, ...
    eur_swap_payments, startDate, dates, discounts);

% Display the European Swaption convergence table and final error
disp(eur_table);
fprintf('Target Analytical Price (Jamshidian):  %.6f\n', eur_analytical_price);
fprintf('Error at 365 steps:                    %.2e\n\n', abs(eur_tree_prices(end) - eur_analytical_price));


% POINT C: BOUNDS VERIFICATION
% Calculate Analytical Bounds
[LB, UB_eur, UB_cap] = compute_bermudan_bounds(a, sigma, K, exercise_dates, ...
    swap_payment_dates, startDate, dates, discounts);

summary_table = verify_swaption_bounds(precision_levels, bermudan_prices, LB, UB_eur, UB_cap);

% Display the results
disp(summary_table);
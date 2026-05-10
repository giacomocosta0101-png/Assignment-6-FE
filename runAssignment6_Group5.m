% runAssignment6
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

%% Case Study 2: Bermudan Swaption Pricing via Hull-White Model

% 1. SETTINGS & MARKET DATA 
settlement_date = datetime(2008, 2, 15);

% Calculate the start date as Settlement + 2 business days
startDate_num = business_date_offset(settlement_date, 'day_offset', 2);
startDate = datetime(startDate_num, 'ConvertFrom', 'datenum');

% Hull-White Model parameters
a = 0.11;       % Mean reversion speed
sigma = 0.008;  % Volatility of the short rate
K = 0.05;       % Fixed strike rate (coupon) of the underlying swap

% Convergence settings: test different number of time steps per year
precision_levels = [1, 4, 12, 52, 365]; 
prices = zeros(length(precision_levels), 1);

% 2. KEY DATES PREPARATION
% Generate Exercise dates: annually from Year 2 to Year 9
y_ex = (2:9)';
exercise_dates = datetime(arrayfun(@(y) business_date_offset(startDate, ...
    'year_offset', y), y_ex), 'ConvertFrom', 'datenum');

% Generate Underlying Swap payment dates: annually from Year 1 to Year 10
y_pay = (1:10)';
swap_payment_dates = datetime(arrayfun(@(y) business_date_offset(startDate,...
    'year_offset', y), y_pay), 'ConvertFrom', 'datenum');

% Final maturity is the last payment date of the swap
maturity_date = swap_payment_dates(end);

% 3. CONVERGENCE LOOP
% Iterate through each precision level to observe price stabilization
for j = 1:length(precision_levels)
    
    % STEP A: GRID CONSTRUCTION
    % Build the refined time grid that includes both uniform steps and exercise dates
    [time_grid_y, grid_dates, dt_vector, idx_exercise] = ...
        build_hw_time_grid(startDate, maturity_date, exercise_dates, precision_levels(j));
    
    % STEP B: TREE GEOMETRY
    % Calculate spatial step (dx) and maximum nodes (l_max) based on average dt
    dt_ref = mean(dt_vector);
    
    % Reference volatility for a single step
    sigma_hat_ref = sigma * sqrt((1 - exp(-2 * a * dt_ref)) / (2 * a));
    dx = sigma_hat_ref * sqrt(3); % Trinomial tree spatial step
    
    % Reference drift adjustment for boundary calculation
    mu_hat_ref = 1 - exp(-a * dt_ref);
    
    % Determine l_max (the level where branching shifts to ensure stability)
    l_max = ceil((1 - sqrt(2/3)) / mu_hat_ref );
    
    % STEP C: PRICING 
    prices(j) = price_bermudan(a, sigma, K, dt_vector, grid_dates, l_max, dx, ...
                                    startDate, dates, discounts, ...
                                    swap_payment_dates, idx_exercise);
    
    % STEP D: OUTPUT RESULTS 
    fprintf('Steps/Year: %3d | Nodes: %4d | Price: %.6f\n', ...
            precision_levels(j), length(grid_dates), prices(j));
end

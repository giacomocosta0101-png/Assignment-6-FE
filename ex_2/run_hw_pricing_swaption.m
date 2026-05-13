function [prices, results_table] = run_hw_pricing_swaption(a, sigma, K, ...
    startDate, exerciseDates, swapPaymentDates, precision_levels, ...
    curve_dates, curve_discounts, instrument_type)
% RUN_HW_PRICING_SWAPTION Manages the Bermudan Swaption pricing workflow.
%
% This function manages the grid construction, tree parameter 
% calculation, and the convergence loop for pricing a Bermudan Swaption 
% using the Hull-White trinomial tree.
%
% INPUTS:
%   a                  : [scalar] Speed of mean reversion in the Hull-White model.
%   sigma              : [scalar] Volatility parameter in the Hull-White model.
%   K                  : [scalar] Fixed strike rate of the swap.
%   startDate          : [datetime/datenum] Valuation date (t=0).
%   maturityDate       : [datetime/datenum] Final maturity of the swap.
%   exerciseDates      : [vector] Possible exercise dates for the swaption.
%   swapPaymentDates   : [vector] Payment dates of the underlying swap.
%   precision_levels   : [vector] Steps per year to test for convergence.
%   curve_dates        : [vector] Dates of the market discount curve.
%   curve_discounts    : [vector] Discount factors from the market curve.
%   instrument_type    : [string] 'Bermudan', 'European', or 'ZCB'.
%
% OUTPUTS:
%   prices             : [vector] Bermudan prices for each precision level.
%   results_table      : [table] Summary of steps, nodes, and prices.

    % Initialization
    n_levels = length(precision_levels);
    prices = zeros(n_levels, 1);
    num_nodes = zeros(n_levels, 1);
    maturityDate = swapPaymentDates(end);

    % CONVERGENCE LOOP
    for j = 1:n_levels
        % 1. GRID CONSTRUCTION
        % Build the refined time grid (snapping to exercise dates)
        [~, grid_dates, dt_vector, idx_exercise] = ...
            build_hw_time_grid(startDate, maturityDate, precision_levels(j), exerciseDates);
        
        num_nodes(j) = length(grid_dates);

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

        % 3. PRICING
        % Select the appropriate pricing engine based on the instrument type
        if strcmpi(instrument_type, 'ZCB')
            % Call the Zero Coupon Bond function (ignores K and exercise dates)
            prices(j) = price_zcb_tree(a, sigma, dt_vector, grid_dates, l_max, dx, ...
                                       startDate, curve_dates, curve_discounts);
                                       
        elseif strcmpi(instrument_type, 'Bermudan') || strcmpi(instrument_type, 'European')
            % Call the standard Swaption pricing function
            prices(j) = price_bermudan(a, sigma, K, dt_vector, grid_dates, l_max, dx, ...
                                       startDate, curve_dates, curve_discounts, ...
                                       swapPaymentDates, idx_exercise);
        else
            error('Invalid instrument type. Use: ''ZCB'', ''Bermudan'', or ''European''.');
        end
    end

    % Create a summary table
    results_table = table(precision_levels(:), num_nodes, prices, ...
        'VariableNames', {'Steps_Per_Year', 'Total_Nodes', 'Bermudan_Price'});
end
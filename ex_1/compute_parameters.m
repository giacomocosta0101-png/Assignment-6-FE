function [parameters,sigma,kappa,eta] = compute_parameters(data, dates, ...
    discounts, strike, spread,volatilities,strikes)
% Pre-computes all the date / discount / forward quantities needed to price
% the certificate swap (Annex 2).
%
% Key choices:
%   - The forward F0 is computed to the COUPON RESET DATE (where S is
%     observed for the digital), not to the payment date. Same horizon will
%     drive the NIG dynamics, so day count is consistent (Act/365).
%   - Coupon accruals use 30/360 between consecutive payment dates, per
%     termsheet ("Payable annually on a 30/360 day basis").
%   - Both adjusted and unadjusted payment dates are stored: BPV needs the
%     unadjusted ones as anchor for the quarterly back-step.

    spot_price = data.reference;
    dividend   = data.dividends;
    
    DC_ACT365 = 3;
    DC_30_360 = 6;
    
    % --- Date pillars -----------------------------------------------------
    % t0 of the swap = valuation + 2bd spot lag (NOT the valuation date)
    valuationDate = datenum(data.date, 'yyyymmdd');
    startDate     = business_date_offset(valuationDate, 'day_offset', 2);
    
    % Adjusted coupon payment dates (modified following)
    firstCouponDate = business_date_offset(startDate, 'year_offset', 1, ...
                          'convention', 'modified_following');
    maturityDate    = business_date_offset(startDate, 'year_offset', 2, ...
                          'convention', 'modified_following');
    
    % Unadjusted payment dates: anchor for BPV's quarterly back-step.
    % Use calyears (= "same day, next year"), NOT years (= 365.2425 days).
    startDate_dt          = datetime(startDate, 'ConvertFrom', 'datenum');
    firstCouponDate_unadj = datenum(startDate_dt + calyears(1));
    maturityDate_unadj    = datenum(startDate_dt + calyears(2));
    
    % Reset date = 2bd before each payment (per termsheet). Spot is observed here.
    couponResetDate = business_date_offset(firstCouponDate, 'day_offset', -2, ...
                          'convention', 'modified_following');
    
    % --- Discount factors -------------------------------------------------
    B_t1    = get_discount_factor_by_zero_rates_linear_interp(startDate, ...
                  firstCouponDate, dates, discounts);
    B_t2    = get_discount_factor_by_zero_rates_linear_interp(startDate, ...
                  maturityDate,    dates, discounts);
    B_reset = get_discount_factor_by_zero_rates_linear_interp(startDate, ...
                  couponResetDate, dates, discounts);
    
    % --- Year fractions ---------------------------------------------------
    % Act/365 for everything tied to the curve and to the model dynamics.
    ttm_reset = yearfrac(startDate, couponResetDate, DC_ACT365);
    
    % 30/360 for the contractual coupon accrual factors. Each spans the
    % period between two CONSECUTIVE PAYMENT dates (not involving the reset).
    accrual_1 = yearfrac(startDate,       firstCouponDate, DC_30_360);
    accrual_2 = yearfrac(firstCouponDate, maturityDate,    DC_30_360);
    
    % --- Forward to the observation date ----------------------------------
    % Continuous-compounded zero rate at the reset horizon
    rate_reset = -log(B_reset) / ttm_reset;
    
    % F0 to the reset date — same horizon and day count as the dynamics
    F0_reset = spot_price * exp((rate_reset - dividend) * ttm_reset);
    
    % Calibration: 
    M = 15;                         % Exponent for the number of points in FFT
    N = 2^M;                        % Total number of points for the FFT grid
    x   = log(F0_reset ./ strikes); % Log-moneyness
    dx = 0.0025;                    % Step size for x-discretization (used in FFT)
    alpha = 0.5;

    p0 = [0.15, 0.5, 10];  % Initial guess p0 = [sigma, kappa, eta]
    lb = [0.01, 0.1, 3];  % Lower bounds
    ub = [1.0, 1.0, 50];   % Upper bounds

    % Solver options: display progress and set convergence tolerance
    options = optimoptions('lsqnonlin', 'Display', 'iter', 'FunctionTolerance', 1e-6);

    % We run the global calibration (weights = 1)
    p_opt = lsqnonlin(@(p) calibration_residuals(p, x, volatilities,...
    ttm_reset, B_reset, F0_reset, rate_reset, alpha, N, dx), p0, lb, ub, options);
    sigma = p_opt(1); kappa = p_opt(2); eta = p_opt(3);

    sigma_black = interp1(strikes,volatilities,strike,'spline');
    
    
    % --- Pack -------------------------------------------------------------
    parameters = struct( ...
        'spot',                  spot_price,            ...
        'dividend',              dividend,              ...
        'strike',                strike,                ...
        'spread',                spread,                ...
        'valuationDate',         valuationDate,         ...
        'startDate',             startDate,             ...
        'firstCouponDate',       firstCouponDate,       ...
        'maturityDate',          maturityDate,          ...
        'firstCouponDate_unadj', firstCouponDate_unadj, ...
        'maturityDate_unadj',    maturityDate_unadj,    ...
        'couponResetDate',       couponResetDate,       ...
        'B_t1',                  B_t1,                  ...
        'B_t2',                  B_t2,                  ...
        'ttm_reset',             ttm_reset,             ...
        'accrual_1',             accrual_1,             ...
        'accrual_2',             accrual_2,             ...
        'F0_reset',              F0_reset,              ...
        'dates_bootstrap',       dates,                 ...
        'discounts_bootstrap',   discounts,             ...
        'sigma_black',           sigma_black,           ...
        'volatilities',          volatilities, ...
        'strikes',               strikes, ...
        'rate',                  rate_reset,...
        'B_reset',                B_reset...
        );
end
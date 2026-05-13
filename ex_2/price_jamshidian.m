function price = price_jamshidian(a, sigma, K, t_exercise, swap_pay_dates, ...
                                  start_date, curve_dates, curve_discounts)
% PRICE_JAMSHIDIAN Exact price of a European Swaption using Jamshidian Decomposition.
%
% This function prices a European Payer Swaption under the Hull-White model.
% It decomposes the option into a portfolio of Put options on Zero Coupon Bonds (ZCB).
%
% INPUTS:
%   a                  : [scalar] Speed of mean reversion in the Hull-White model.
%   sigma              : [scalar] Volatility parameter in the Hull-White model.
%   K                  : [scalar] Fixed strike rate of the underlying swap 
%   t_exercise         : [scalar] Option expiry date. The date when you decide to exercise 
%                                 the swaption (also the Start Date of the underlying swap).
%   swap_pay_dates     : [vector] Dates of the swap cash flows (coupons and principal).
%   start_date         : [scalar] Today's date (t=0) for discounting.
%   curve_dates        : [vector] Dates of the market discount curve.
%   curve_discounts    : [vector] Discount factors from the market curve.
%
% OUTPUT:
%   price              : [scalar] European Swaption price at t=0.

    % 1. MARKET DATA & YEAR FRACTIONS
    % Market ZCB from t=0 to exercise and to each payment date
    B0_T_ex = get_discount_factor_by_zero_rates_linear_interp(start_date, ...
                t_exercise, curve_dates, curve_discounts);
    B0_T_pay = get_discount_factor_by_zero_rates_linear_interp(start_date, ...
                swap_pay_dates, curve_dates, curve_discounts);
    
    % Time in years from t=0 (Act/365 convention for HW formulas)
    T_ex_y = yearfrac(start_date, t_exercise, 3);
    T_j_y = yearfrac(start_date, swap_pay_dates, 3);
    
    % Fixed leg accruals (30/360 convention)
    tau_pay = [t_exercise; swap_pay_dates(:)];
    accruals = yearfrac(tau_pay(1:end-1), tau_pay(2:end), 6); 
    
    % Hull-White B coefficient: B(t_ex, T_j) = (1 - exp(-a*(T_j - t_ex)))/a
    B_ex_pay = (1 - exp(-a * (T_j_y - T_ex_y))) / a;
    
    % 2. FIND x*
    % x is the stochastic component of the short rate: r(t) = x(t) + phi(t)
    % We search for x* such that the Value of the Swap at t_exercise is ZERO.
    
    % Cash flow weights (Coupons c_j and Final Principal at the last date)
    c_j = K * accruals;
    c_j(end) = c_j(end) + 1;
    
    % Forward Discount Factor: B(0, T_j) / B(0, t_ex)
    Fwd_ZCB = B0_T_pay ./ B0_T_ex;
    
    % HW Convexity adjustment V(t_ex, T_j)
    V_t_T = (sigma^2 / (2*a^3)) * (1 - exp(-a*(T_j_y - T_ex_y))).^2 * (1 - exp(-2*a*T_ex_y));
    
    % Objective function: Sum_{j} c_j * P(t_ex, T_j, x) = 1
    % This represents the condition for the swaption to be At-The-Money at exercise.
    obj_fun = @(x) sum( c_j .* Fwd_ZCB .* exp(-B_ex_pay * x - 0.5 * V_t_T) ) - 1;
    
    % x_star is the value of the stochastic state that makes the swap par
    x_star = fzero(obj_fun, 0.0);
    
    % 3. PRICE INDIVIDUAL ZCB OPTIONS
    % Jamshidian Strike K_j for each ZCB: K_j = B(t_ex, T_j, x*)
    K_j = Fwd_ZCB .* exp(-B_ex_pay * x_star - 0.5 * V_t_T);
    
    % Integrated volatility sigma_p for the ZCB option
    sigma_p = (sigma/a) * (1 - exp(-a*(T_j_y - T_ex_y))) * sqrt((1 - exp(-2*a*T_ex_y))/(2*a));
    
    % Black formula for Put option on ZCB
    d1 = (log(Fwd_ZCB ./ K_j) + 0.5 * sigma_p.^2) ./ sigma_p;
    d2 = d1 - sigma_p;
    
    % Put Price = P(0, t_ex) * [ K_j * N(-d2) - Fwd_ZCB * N(-d1) ]
    zcb_options = B0_T_ex * (K_j .* normcdf(-d2) - Fwd_ZCB .* normcdf(-d1));
    
    % The Swaption price is the sum of ZCB options weighted by cash flows c_j
    price = sum(c_j .* zcb_options);
end
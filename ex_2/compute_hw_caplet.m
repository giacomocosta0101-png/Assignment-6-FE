function caplet_price = compute_hw_caplet(a, sigma, K, T_reset, T_payment, ...
    start_date, curve_dates, curve_discounts)
% COMPUTE_HW_CAPLET Analytical price of a single caplet under Hull-White.
%
% Prices a caplet for period [T_reset, T_payment] with strike K using the
% standard equivalence between a caplet and a put option on a Zero Coupon
% Bond:
%   Caplet = (1 + K*alpha) * Put_on_ZCB(K* = 1/(1+K*alpha), T_reset, T_payment)
% where Put_on_ZCB is the Black formula on the forward ZCB B(T_reset, T_payment)
% with maturity T_reset and integrated Hull-White volatility sigma_p.
%
% Used as a building block for the tight Bermudan upper bound:
%   V_Berm <= sum_i Caplet_i
% one caplet per future period reachable by the Bermudan exercise schedule.
%
% INPUTS:
%   a, sigma           : [scalar] Hull-White model parameters.
%   K                  : [scalar] Strike rate (same as Bermudan strike).
%   T_reset            : [scalar] Libor reset date (also caplet maturity), datenum.
%   T_payment          : [scalar] Caplet payment date, datenum.
%   start_date         : [scalar] Valuation date t=0, datenum.
%   curve_dates        : [vector] Market discount curve dates.
%   curve_discounts    : [vector] Market discount factors.
%
% OUTPUT:
%   caplet_price       : [scalar] Caplet price at t=0.

    % Year fractions (Act/365 for HW dynamics)
    T_reset_y   = yearfrac(start_date, T_reset,   3);
    T_pay_y     = yearfrac(start_date, T_payment, 3);

    % Accrual factor on the underlying Libor period (30/360, as for fixed leg)
    alpha = yearfrac(T_reset, T_payment, 6);

    % Market discount factors P(0, T_reset) and P(0, T_payment)
    B0_reset = get_discount_factor_by_zero_rates_linear_interp(start_date, ...
                T_reset,   curve_dates, curve_discounts);
    B0_pay   = get_discount_factor_by_zero_rates_linear_interp(start_date, ...
                T_payment, curve_dates, curve_discounts);

    % Forward ZCB at T_reset for delivery at T_payment
    F_ZCB = B0_pay / B0_reset;

    % Equivalent strike on the ZCB: K* = 1 / (1 + K * alpha)
    K_star = 1 / (1 + K * alpha);

    % Hull-White integrated volatility for the bond option
    % sigma_p = (sigma/a) * (1 - exp(-a*tau)) * sqrt((1 - exp(-2*a*T_reset)) / (2*a))
    tau     = T_pay_y - T_reset_y;
    sigma_p = (sigma / a) * (1 - exp(-a * tau)) * ...
              sqrt((1 - exp(-2 * a * T_reset_y)) / (2 * a));

    % Black formula for a Put on the ZCB
    d1 = (log(F_ZCB / K_star) + 0.5 * sigma_p^2) / sigma_p;
    d2 = d1 - sigma_p;
    put_zcb = B0_reset * (K_star * normcdf(-d2) - F_ZCB * normcdf(-d1));

    % Caplet = leverage factor * put on ZCB
    caplet_price = (1 + K * alpha) * put_zcb;
end

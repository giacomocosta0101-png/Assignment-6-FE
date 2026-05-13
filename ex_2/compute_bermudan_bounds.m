function [LB, UB_eur, UB_cap] = compute_bermudan_bounds(a, sigma, K, ...
    exercise_dates, swap_pay_dates, start_date, curve_dates, curve_discounts)
% COMPUTE_BERMUDAN_BOUNDS Analytical bounds for a Bermudan Swaption.
%
% Returns three bounds that the Bermudan price V_Berm must satisfy:
%
%   LB  = max_i V_Eur_i (residual swap)                  -- lower bound (tight)
%   UB_eur = sum_i V_Eur_i (residual swap)               -- upper bound (loose)
%   UB_cap = sum_i Caplet_i (single-period [T_i, T_i+1]) -- upper bound (tight)
%
% THEORY:
%
% LB: A Bermudan can always replicate any single-exercise European policy,
%     so V_Berm >= max_i V_Eur_i.
%
% UB_eur: Decomposing the Bermudan as a strip of independent exercise rights
%     on the full residual swap, each right is at most as valuable as the
%     corresponding European on the residual swap; summing gives an upper
%     bound. Loose because cash flows after T_i are double-counted across i.
%
% UB_cap: Decomposing the swap payoff per cash flow and applying (sum x)^+ <=
%     sum (x^+), one obtains
%        V_Berm <= sum_j Caplet[T_{j-1}, T_j]
%     i.e. a strip of single-period caplets, one per future swap period
%     reachable by some Bermudan exercise. Each cash flow is counted exactly
%     once, so UB_cap is typically 4-5x tighter than UB_eur.
%     Each caplet is priced analytically under Hull-White via the put-on-ZCB
%     equivalence.
%
% INPUTS:
%   a, sigma           : Hull-White parameters.
%   K                  : Strike of the underlying swap.
%   exercise_dates     : Bermudan exercise dates (datenum).
%   swap_pay_dates     : Full swap payment schedule (datenum).
%   start_date         : Valuation date t=0 (datenum).
%   curve_dates        : Market discount curve dates.
%   curve_discounts    : Market discount factors.
%
% OUTPUTS:
%   LB      : Lower bound (max of European swaptions).
%   UB_eur  : Loose upper bound (sum of European swaptions).
%   UB_cap  : Tight upper bound (sum of single-period caplets).

    n_ex          = length(exercise_dates);
    eur_prices    = zeros(n_ex, 1);
    caplet_prices = zeros(n_ex, 1);

    for i = 1:n_ex
        t_ex = exercise_dates(i);

        % European on the full residual swap (for LB and UB_eur) 
        % When exercising at t_ex, only payments strictly after t_ex are received.
        relevant_payments = swap_pay_dates(swap_pay_dates > t_ex);
        eur_prices(i) = price_jamshidian(a, sigma, K, t_ex, relevant_payments, ...
                                          start_date, curve_dates, curve_discounts);

        % Single-period caplet [t_ex, next swap payment] (for UB_cap)
        % Each Bermudan exercise locks in at most one immediate-period cash flow;
        % the rest of the residual swap is bounded by later caplets.
        next_payment       = min(swap_pay_dates(swap_pay_dates > t_ex));
        caplet_prices(i)   = compute_hw_caplet(a, sigma, K, t_ex, next_payment, ...
                                                start_date, curve_dates, curve_discounts);
    end

    % Bounds
    LB     = max(eur_prices);
    UB_eur = sum(eur_prices);
    UB_cap = sum(caplet_prices);
end

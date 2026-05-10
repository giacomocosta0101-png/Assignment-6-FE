function payoff = calc_swap_payoff(x_nodes, t_exercise, K, a, sigma, ...
        T_payments, B_t0_t_exercise, B_t0_t_payments, start_date)
% CALC_SWAP_PAYOFF Computes the payer swap intrinsic value across all grid nodes.
%
% This function calculates the value of the underlying swap at an exercise 
% date for each node of the trinomial tree using the Hull-White closed-form 
% formula for Zero Coupon Bonds.
%
% INPUTS:
%   x_nodes           : [vector] Spatial grid nodes (x = r - alpha).
%   t_exercise        : [scalar] Current exercise date (datenum).
%   K                 : [scalar] Strike rate (fixed coupon).
%   a                 : [scalar] Speed of mean reversion in the Hull-White model.
%   sigma             : [scalar] Volatility parameter in the Hull-White model.
%   T_payments        : [vector] Future swap payment dates.
%   B_t0_t_exercise   : [scalar] Market discount factor P(0, t_exercise).
%   B_t0_t_payments   : [vector] Market discount factors P(0, T_payments).
%   start_date        : [scalar] Valuation date (t=0).
%
% OUTPUT:
%   payoff            : [vector] Max(0, Swap Value) for each node.

    % 1. DATA PREPARATION
    x_nodes = x_nodes(:);             % [N_nodes x 1]
    T_payments = T_payments(:);       % [N_payments x 1]
    B_t0_t_payments = B_t0_t_payments(:); % [N_payments x 1]
    
    % Time from start to exercise (for Hull-White variance component)
    yfrac_t0_t_ex = yearfrac(start_date, t_exercise, 3);
    V_t = (sigma^2 / (2 * a)) * (1 - exp(-2 * a * yfrac_t0_t_ex));
    
    % 2. ACCRUAL FACTORS (Fixed Leg - 30/360 convention)
    all_dates = [t_exercise; T_payments];
    accrual_j = yearfrac(all_dates(1:end-1), all_dates(2:end), 6); 
    
    % 3. HULL-WHITE BOND COMPONENTS (B_tau and Variance Integral)
    tau_j = yearfrac(t_exercise, T_payments, 3);
    B_tau = (1 - exp(-a * tau_j)) ./ a;  % [N_payments x 1]
    var_adj = 0.5 * V_t .* (B_tau.^2);    % [N_payments x 1]
    
    % 4. ZERO COUPON BOND MATRIX (ZCB)
    % We compute B(t_exercise, T_j) for each node and each payment date.
    % Dimensions: [N_nodes x N_payments]
    % Formula: P(0,T)/P(0,t) * exp( -B_tau*x_nodes - 0.5*B_tau^2*V_t )
    
    % Forward market discount factors
    forward_ZCB = (B_t0_t_payments ./ B_t0_t_exercise)'; % [1 x N_payments]
    
    % Exponent calculation using broadcasting (x_nodes * B_tau')
    % Result is [N_nodes x N_payments]
    exponent = -(x_nodes * B_tau') - var_adj';
    
    % Full ZCB matrix across all nodes and all future payment dates
    ZCB_matrix = forward_ZCB .* exp(exponent);
    
    % 5. SWAP VALUE CALCULATION
    % Fixed Leg: sum of (accrual * K * ZCB) for each payment
    % Vector-Matrix product results in [N_nodes x 1]
    Fixed_Leg = ZCB_matrix * (accrual_j .* K);
    
    % Floating Leg: 1 - ZCB_maturity (Principle of par swap)
    V_float = 1 - ZCB_matrix(:, end);
    
    % 6. PAYOFF (Payer Swaption)
    % Value = Max(0, Floating Leg - Fixed Leg)
    payoff = max(0, V_float - Fixed_Leg);
    
end
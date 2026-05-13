function payoff = calc_swap_payoff(x_nodes, t_exercise, K, a, sigma, ...
        T_payments, B_t0_t_exercise, B_t0_t_payments, start_date)
% CALC_SWAP_PAYOFF Computes the payer swap intrinsic value across all grid nodes.
%
% This function calculates the value of the underlying swap at an exercise 
% date for each node of the trinomial tree using the Hull-White closed-form 
% formula for Zero Coupon Bonds (via compute_hw_zcb).

    % 1. DATA PREPARATION
    T_payments = T_payments(:);       % Ensure column vector
    
    % 2. ACCRUAL FACTORS (Fixed Leg - 30/360 convention)
    all_dates = [t_exercise; T_payments];
    accrual_j = yearfrac(all_dates(1:end-1), all_dates(2:end), 6); 
    
    % 3. CALL ANALYTICAL ZCB FUNCTION (The Engine)
    % Qui usiamo la nuova funzione robusta che abbiamo appena sistemato!
    % Restituirà la matrice esatta [N_nodes x N_payments] senza NaN.
    ZCB_matrix = compute_hw_zcb(x_nodes, t_exercise, T_payments, a, sigma, ...
                                B_t0_t_exercise, B_t0_t_payments, start_date);
    
    % 4. SWAP VALUE CALCULATION
    % Fixed Leg: sum of (accrual * K * ZCB) for each payment
    % Prodotto Vettore-Matrice: [N_nodes x N_payments] * [N_payments x 1] -> [N_nodes x 1]
    Fixed_Leg = ZCB_matrix * (accrual_j .* K);
    
    % Floating Leg: 1 - ZCB_maturity (Principle of par swap)
    V_float = 1 - ZCB_matrix(:, end);
    
    % 5. PAYOFF (Payer Swaption)
    % Value = Max(0, Floating Leg - Fixed Leg)
    payoff = max(0, V_float - Fixed_Leg);
    
end
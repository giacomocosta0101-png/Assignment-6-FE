function price = price_zcb_tree(a, sigma, dt_vector, grid_dates, l_max, dx, ...
                                start_date, curve_dates, curve_discounts)
% PRICE_ZCB_TREE Prices a Zero Coupon Bond using the Hull-White Trinomial Tree.
%
% This function is used to test the "perfect calibration" of the tree.
% The output price should exactly match the deterministic market discount 
% factor P(0, T_maturity) extracted from the initial yield curve.

    % 1. DATA TYPE UNIFICATION
    grid_dates = datenum(grid_dates);
    start_date = datenum(start_date);
    curve_dates = datenum(curve_dates);

    % 2. GRID INITIALIZATION
    l = (-l_max : l_max)';          
    N_nodes = length(l);
    N_steps = length(dt_vector);
    x_grid = l * dx;                
    
    % TERMINAL CONDITION: A ZCB pays exactly 1 at maturity in all states
    V = ones(N_nodes, 1); 

    % 3. BACKWARD INDUCTION LOOP (Pure linear discounting)
    for i = N_steps : -1 : 1
        dt = dt_vector(i);
        t_curr = grid_dates(i);
        t_next = grid_dates(i+1);
        
        % Tree parameters & Probabilities
        mu_hat = 1 - exp(-a * dt);
        M = l * mu_hat;
        
        pu = zeros(N_nodes, 1); pm = zeros(N_nodes, 1); pd = zeros(N_nodes, 1);
        idx_A = 2 : N_nodes - 1;    
        
        pu(idx_A) = 1/6 + 0.5*(M(idx_A).^2 - M(idx_A));
        pm(idx_A) = 2/3 - M(idx_A).^2;
        pd(idx_A) = 1/6 + 0.5*(M(idx_A).^2 + M(idx_A));
        
        pu(1) = 1/6 + 0.5*(M(1)^2 + M(1)); 
        pm(1) = -1/3 - M(1)^2 - 2*M(1); 
        pd(1) = 7/6 + 0.5*(M(1)^2 + 3*M(1));
        
        pu(end) = 7/6 + 0.5*(M(end)^2 - 3*M(end)); 
        pm(end) = -1/3 - M(end)^2 + 2*M(end); 
        pd(end) = 1/6 + 0.5*(M(end)^2 - M(end));

        % 4. DETERMINISTIC FORWARD SHIFT (The crucial step we discussed!)
        B0_t = get_discount_factor_by_zero_rates_linear_interp(start_date,...
                t_curr, curve_dates, curve_discounts);
        B0_T = get_discount_factor_by_zero_rates_linear_interp(start_date,...
                t_next, curve_dates, curve_discounts);

        B_HW_analytical = compute_hw_zcb(x_grid, t_curr, t_next, a, sigma, ...
                B0_t, B0_T, start_date);
               
        % 5. STOCHASTIC DISCOUNT
        [Du, Dm, Dd] = compute_stochastic_discount(x_grid, B_HW_analytical, dx, mu_hat, a, sigma, dt);
                       
        % 6. CONTINUATION VALUE
        V_cont = zeros(N_nodes, 1);
        V_cont(idx_A) = Du(idx_A).*pu(idx_A).*V(idx_A+1) + ...
                        Dm(idx_A).*pm(idx_A).*V(idx_A)   + ...
                        Dd(idx_A).*pd(idx_A).*V(idx_A-1);
        
        V_cont(1)     = Dd(1)*pd(1)*V(1) + Dm(1)*pm(1)*V(2) + Du(1)*pu(1)*V(3);
        V_cont(end)   = Du(end)*pu(end)*V(end) + Dm(end)*pm(end)*V(end-1) + Dd(end)*pd(end)*V(end-2);
        
        % NO EARLY EXERCISE (It's a simple bond, so V is just V_cont)
        V = V_cont;
    end

    price = V(l == 0); 
end
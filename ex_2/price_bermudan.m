function price = price_bermudan(a, sigma, K, dt_vector, grid_dates, l_max, dx, ...
    start_date, curve_dates, curve_discounts, ...
    swap_payment_dates, idx_exercise)
% PRICE_BERMUDAN Prices a Bermudan Swaption using a vectorized Hull-White Trinomial Tree.
%
% This function performs backward induction on a trinomial tree with variable 
% time steps, incorporating Hull-White stochastic discounting and early exercise.
%
% INPUTS:
%   a                  : [scalar] Speed of mean reversion in the Hull-White model.
%   sigma              : [scalar] Volatility parameter in the Hull-White model.
%   K                  : [scalar] Fixed strike rate of the underlying swap.
%   dt_vector          : [vector] Time intervals between grid nodes [N_steps x 1].
%   grid_dates         : [vector] All dates in the tree grid [N_steps+1 x 1] (datetime/datenum).
%   l_max              : [scalar] Maximum number of spatial nodes in each direction.
%   dx                 : [scalar] Spatial step size (Delta x).
%   start_date         : [scalar] Valuation/Start date t=0 (datetime or datenum).
%   curve_dates        : [vector] Dates of the market discount curve [N_curve x 1].
%   curve_discounts    : [vector] Discount factors from the market curve [N_curve x 1].
%   swap_payment_dates : [vector] Payment schedule of the underlying fixed leg [N_payments x 1].
%   idx_exercise       : [vector] Indices of exercise dates within grid_dates [N_ex x 1].
%
% OUTPUT:
%   price              : [scalar] Present value of the Bermudan swaption at t=0.

    % 1. DATA TYPE UNIFICATION
    % Convert all date inputs to datenum for high-speed numerical interpolation
    grid_dates = datenum(grid_dates);
    swap_payment_dates = datenum(swap_payment_dates);
    start_date = datenum(start_date);
    curve_dates = datenum(curve_dates);

    % 2. GRID INITIALIZATION
    l = (-l_max : l_max)';          % Spatial node indices (column vector)
    N_nodes = length(l);
    N_steps = length(dt_vector);
    x_grid = l * dx;                % Pre-compute spatial grid values (x = r - alpha)
    
    % Option value at maturity is 0 (terminal condition)
    V = zeros(N_nodes, 1); 

    % 3. BACKWARD INDUCTION LOOP
    for i = N_steps : -1 : 1
        dt = dt_vector(i);
        t_curr = grid_dates(i);
        t_next = grid_dates(i+1);
        
        % Tree parameters
        mu_hat = 1 - exp(-a * dt);
        M = l * mu_hat;

        % Initialization of Trinomial Tree Probabilities 
        pu = zeros(N_nodes, 1); 
        pm = zeros(N_nodes, 1); 
        pd = zeros(N_nodes, 1);
        idx_A = 2 : N_nodes - 1;    % Indices for standard central branching
        
        % CASE A: Standard branching (Central nodes)
        pu(idx_A) = 1/6 + 0.5*(M(idx_A).^2 - M(idx_A));
        pm(idx_A) = 2/3 - M(idx_A).^2;
        pd(idx_A) = 1/6 + 0.5*(M(idx_A).^2 + M(idx_A));
        
        % CASE B: Bottom boundary (Node 1) - Upward branching
        pu(1) = 1/6 + 0.5*(M(1)^2 + M(1)); 
        pm(1) = -1/3 - M(1)^2 - 2*M(1); 
        pd(1) = 7/6 + 0.5*(M(1)^2 + 3*M(1));
        
        % CASE C: Top boundary (Last node) - Downward branching
        pu(end) = 7/6 + 0.5*(M(end)^2 - 3*M(end)); 
        pm(end) = -1/3 - M(end)^2 + 2*M(end); 
        pd(end) = 1/6 + 0.5*(M(end)^2 - M(end));

        % 4. MARKET DATA & STOCHASTIC DISCOUNTING
        % Get deterministic discount factors between current and next node
        B0_t = get_discount_factor_by_zero_rates_linear_interp(start_date,...
                t_curr, curve_dates, curve_discounts);
        B0_T = get_discount_factor_by_zero_rates_linear_interp(start_date,...
                t_next, curve_dates, curve_discounts);
    
        B_HW_analytical = compute_hw_zcb(x_grid, t_curr, t_next, a, sigma, ...
                B0_t, B0_T, start_date);
               
       % 5. STOCHASTIC DISCOUNT
        [Du, Dm, Dd] = compute_stochastic_discount(x_grid, B_HW_analytical, dx, mu_hat, a, sigma, dt);
        
       % CONTINUATION VALUE CALCULATION
        V_cont = zeros(N_nodes, 1);
        
        % Apply branching logic to continuation value (Vectorized)
        V_cont(idx_A) = Du(idx_A).*pu(idx_A).*V(idx_A+1) + ...
                        Dm(idx_A).*pm(idx_A).*V(idx_A)   + ...
                        Dd(idx_A).*pd(idx_A).*V(idx_A-1);
        
        % Edge cases: Boundary nodes use specific upward/downward shifts
        V_cont(1)     = Dd(1)*pd(1)*V(1) + Dm(1)*pm(1)*V(2) + Du(1)*pu(1)*V(3);

        V_cont(end)   = Du(end)*pu(end)*V(end) + Dm(end)*pm(end)*V(end-1) +...
                        Dd(end)*pd(end)*V(end-2);
        
        % 6. EARLY EXERCISE LOGIC (Bermudan feature)
        if any(idx_exercise == i)
            % Identify future swap payments strictly after current exercise date
            T_pay_future = swap_payment_dates(swap_payment_dates > t_curr);
            
            % Market discount factors for future payments
            B0_T_pay = get_discount_factor_by_zero_rates_linear_interp(...
                    start_date, T_pay_future, curve_dates, curve_discounts);
            
            % Compute swap payoff across all nodes (intrinsic value)
            Payoff = calc_swap_payoff(x_grid, t_curr, K, a, sigma, ...
                    T_pay_future, B0_t, B0_T_pay, start_date);
            
            % Optimal decision: Max of continuation or immediate exercise
            V = max(V_cont, Payoff(:));
        else
            % European step: only continuation value
            V = V_cont;
        end
    end

    % Final result is the value at the center of the tree (spatial index l = 0)
    price = V(l == 0); 
end
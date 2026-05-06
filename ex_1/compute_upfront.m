function upfront = compute_upfront(parameters, sigma, kappa, eta, ...
    N_sim,alpha,first_coupon,second_coupon,spread,notional)
% PRICE_MC_NMV Prices Call options using Monte Carlo simulation 
%
% INPUTS:
%   x_grid         - [Vector] Log-moneyness grid: x = ln(F0 / K)
%   sigma          - [Scalar] Volatility parameter 
%   kappa          - [Scalar] Vol-of-vol parameter 
%   eta            - [Scalar] Skewness/Drift adjustment parameter
%   t              - [Scalar] Time to maturity in years
%   N_sim          - [Scalar] Number of Monte Carlo paths to simulate
%   B              - [Scalar] Discount factor for maturity t
%   F0             - [Scalar] Current Forward price
%   alpha          - [Scalar] Tail index (0.5 for NIG, <0.5 for Generalized IG)
%
% OUTPUT:
%   price_MC_NMV   - [Vector] Call prices evaluated at x_grid.
    
    CONVENTION_30_360 = 6;
    F0_EurStoxx_1st_coupon_date = parameters.F0_t1;
    coupon_reset_date = business_date_offset(parameters.firstCouponDate,'day_offset',...
        -2,'convention', 'modified_following');
    accrual_1st_coupon = yearfrac(parameters.startDate,coupon_reset_date,CONVENTION_30_360);
    accrual_2nd_coupon = yearfrac(coupon_reset_date,parameters.maturityDate,CONVENTION_30_360);
    B_1st_coupon_date = parameters.B_t1;
    B_maturity_date = parameters.B_t2;
    start_date = parameters.startDate;
    t = accrual_1st_coupon;
    strike = parameters.strike;
    % Fix seed for reproducibility
    rng(42)

    %% 1. Simulation of the Subordinator G
    % We assume G follows an Inverse Gaussian distribution. 
    % We set the mean 1 (mu_ig = 1) 
    mu_ig = 1; 
    lambda_ig = t / kappa; % lambda_ig is such that the subordinator G has E[G] = 1 And Var(G) = kappa/t 
    % In MATLAB's 'InverseGaussian' parametrization, where Mean = mu and Variance = mu^3/lambda,
    % setting mu = 1 and lambda = t/kappa ensures these statistical properties are met.

    G = random('InverseGaussian', mu_ig, lambda_ig, [N_sim, 1]);
    
    %% 2. Simulation of Standard Normal Variable
    % 'g' represents the independent Gaussian shocks in the mixture model.
    g = randn(N_sim, 1);
    
    %% 3. Martingale Correction (Log-Laplace Exponent)
    % This term (ln_L_eta) is the convexity adjustment. It ensures that 
    % E[exp(ft)] = 1, preserving the Martingale property of the discounted price.
    % The formula incorporates 'alpha' to handle different tail behaviors.
    ln_L_eta = (t / kappa) * ((1 - alpha) / alpha) * ...
               (1 - (1 + (kappa * eta * sigma^2) / (1 - alpha))^alpha);
    
    %% 4. Log-Return Dynamics (ft)
    % The log-return is modeled as a subordinated Brownian motion with drift:
    % ft = diffusion_term - drift_term - martingale_adjustment
    ft = sqrt(t) * sigma .* sqrt(G) .* g - (0.5 + eta) * t * sigma^2 .* G - ln_L_eta;
    
    %% 6 Terminal Price Normalization
    % Calculate the ratio S_T / F0.
    FT_EuroStoxx_1st_coupon_date= F0_EurStoxx_1st_coupon_date*exp(ft); 

    trigger_condition = (FT_EuroStoxx_1st_coupon_date<strike);

    price_bond_sim = notional*(B_1st_coupon_date*trigger_condition*accrual_1st_coupon*...
            first_coupon+B_maturity_date*(1-trigger_condition)*...
            accrual_2nd_coupon*second_coupon);
    price_bond = mean(price_bond_sim,1);
    

    % 4. PV of Party A Leg (Floating Euribor 3M + Spread)
    % PV = Notional * (Floating_Part + Spread_Part)
    first_coupon_date_unadj = start_date +datenum(years(1));
    maturity_date_undadj = start_date +datenum(years(2));% unadjusted maturity
    BPV_early_redemption = BasisPointValueFloating(start_date,... 
        first_coupon_date_unadj, parameters.dates_bootstrap,...
        parameters.discounts_bootstrap);
    BPV_maturity = BasisPointValueFloating(start_date,... 
        maturity_date_undadj, parameters.dates_bootstrap,...
        parameters.discounts_bootstrap);
    
    pv_A_sim = notional*(trigger_condition*((1-B_1st_coupon_date)+spread*BPV_early_redemption)+...
        (1-trigger_condition)*((1-B_maturity_date)+spread*BPV_maturity));
    pv_A = mean(pv_A_sim,1);

    upfront = (pv_A-price_bond)/notional;


end
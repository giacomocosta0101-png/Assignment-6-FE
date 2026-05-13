function upfront = compute_upfront(parameters, sigma, kappa, eta, ...
    N_sim, alpha, first_coupon, second_coupon, spread, notional, ...
    method)
% COMPUTE_UPFRONT
% Mid-market upfront X% for the certificate swap.
%
% Methods:
%   'NIG'        : NIG Monte Carlo
%   'BlackSmile' : Black digital using smile correction through digital_price_smile
%   'BlackFlat'  : Black digital using flat volatility through price_digital_Black

    if nargin < 11 || isempty(method)
        method = 'NIG';
    end

    method = lower(method);

    DC_30_360 = 6;

    % Coupon accruals
    accrual_1 = yearfrac(parameters.startDate, ...
                         parameters.firstCouponDate, DC_30_360);

    accrual_2 = yearfrac(parameters.firstCouponDate, ...
                         parameters.maturityDate, DC_30_360);

    B_t1 = parameters.B_t1;
    B_t2 = parameters.B_t2;

    % Floating leg BPVs
    sd_dt = datetime(parameters.startDate, 'ConvertFrom', 'datenum');

    first_coupon_unadj = parameters.firstCouponDate_unadj;
    maturity_unadj     = parameters.maturityDate_unadj;

    BPV_y1 = BasisPointValueFloating(parameters.startDate, first_coupon_unadj, ...
        parameters.dates_bootstrap, parameters.discounts_bootstrap);

    BPV_y2 = BasisPointValueFloating(parameters.startDate, maturity_unadj, ...
        parameters.dates_bootstrap, parameters.discounts_bootstrap);

    % Trigger probability
    switch method

        case 'nig'

            t  = parameters.ttm_reset;
            F0 = parameters.F0_reset;
            K  = parameters.strike;

            rng(42);

            G = random('InverseGaussian', 1, t/kappa, [N_sim, 1]);
            g = randn(N_sim, 1);

            ln_L_eta = (t/kappa) * ((1 - alpha)/alpha) * ...
                (1 - (1 + (kappa*eta*sigma^2)/(1 - alpha))^alpha);

            ft = sqrt(t)*sigma.*sqrt(G).*g ...
                 - (0.5 + eta)*t*sigma^2.*G ...
                 - ln_L_eta;

            F_t1 = F0 .* exp(ft);

            prob_trigger = mean(F_t1 < K);

        case 'blacksmile'

            % digital_price_smile returns the discounted digital CALL price:
            % price = B_reset * P(S_reset > K)
            dg_price_call = digital_price_smile(parameters);

            prob_no_trigger = dg_price_call / parameters.B_reset;

            prob_trigger = 1 - prob_no_trigger;

        case 'blackflat'

            % price_digital_Black returns discounted put/call digital prices
            price_digital_put= price_digital_Black(parameters, parameters.sigma_black);

            % Trigger is S < K, so use the put digital.
            % price_digital_put = B_reset * P(S_reset < K)
            prob_trigger = price_digital_put / parameters.B_reset;

        otherwise

            error("Unknown method. Use 'NIG', 'BlackSmile', or 'BlackFlat'.");

    end

    prob_no_trigger = 1 - prob_trigger;

    % PV of coupon leg
    pv_coupon = notional * ( ...
        B_t1 * prob_trigger    * accrual_1 * first_coupon + ...
        B_t2 * prob_no_trigger * accrual_2 * second_coupon );

    % PV of Party A floating + spread leg
    pv_A = notional * ( ...
        prob_trigger    * ((1 - B_t1) + spread * BPV_y1) + ...
        prob_no_trigger * ((1 - B_t2) + spread * BPV_y2) );

    upfront = (pv_A - pv_coupon) / notional;

end
function upfront = compute_upfront_3y(parameters, sigma, kappa, eta, ...
    N_sim, alpha, coupon_high, coupon_low, spread, notional)
% Mid-market upfront for the 3y version of the certificate (points c/d).
%
% Payoff: same digital coupon (6%) at year 1 and year 2; early redemption
% fires the first time S falls below K. If neither year triggers, a final
% 2% coupon is paid at year 3.
%
% Point d uses deterministic rates from the curve plus a constant dividend
% yield, with the NIG parameters from point a. The spot under Q reads
%     S_t = S_0 * exp( integrated_drift(t) + f_t )
% so we chain two independent Levy increments to get the spot at the year 1
% and year 2 reset dates.

    rng(42);

    DC_ACT365 = 3;
    DC_30_360 = 6;

    S_0       = parameters.spot;
    div       = parameters.dividend;
    K         = parameters.strike;
    d_curve   = parameters.dates_bootstrap;
    B_curve   = parameters.discounts_bootstrap;
    startDate = parameters.startDate;

    % Payment and reset dates. Year 1 and 2 are reused from `parameters`;
    % only the year 3 dates are new in this case.
    t1_pay   = parameters.firstCouponDate;
    t2_pay   = parameters.maturityDate;
    t3_pay   = business_date_offset(startDate, 'year_offset', 3, ...
                   'convention', 'modified_following');
    t1_reset = parameters.couponResetDate;
    t2_reset = business_date_offset(t2_pay, 'day_offset', -2, ...
                   'convention', 'modified_following');

    % Discount factors at every relevant date.
    B_t1_reset = get_discount_factor_by_zero_rates_linear_interp(startDate, ...
                     t1_reset, d_curve, B_curve);
    B_t2_reset = get_discount_factor_by_zero_rates_linear_interp(startDate, ...
                     t2_reset, d_curve, B_curve);
    B_t1_pay = parameters.B_t1;
    B_t2_pay = parameters.B_t2;
    B_t3_pay = get_discount_factor_by_zero_rates_linear_interp(startDate, ...
                   t3_pay, d_curve, B_curve);

    % Times that drive the dynamics: Act/365, same convention as the curve.
    ttm_1 = parameters.ttm_reset;
    ttm_2 = yearfrac(startDate, t2_reset, DC_ACT365);
    dt_1  = ttm_1;
    dt_2  = ttm_2 - ttm_1;

    % Coupon accruals between consecutive payment dates (30/360 from termsheet).
    accr_1 = parameters.accrual_1;
    accr_2 = parameters.accrual_2;
    accr_3 = yearfrac(t2_pay, t3_pay, DC_30_360);

    % Deterministic part of the spot growth between observations.
    % Rates push S up via 1/B, dividends pull S down via exp(-d * dt).
    grow_1 = exp(-div * dt_1) / B_t1_reset;
    grow_2 = exp(-div * dt_2) * B_t1_reset / B_t2_reset;

    % Two independent NIG increments. Independence is what keeps the
    % multi-period structure consistent with a Levy process; without it
    % we would lose the joint distribution of (S_t1, S_t2).
    Df_1 = nig_increment(sigma, kappa, eta, alpha, dt_1, N_sim);
    Df_2 = nig_increment(sigma, kappa, eta, alpha, dt_2, N_sim);

    % Chain the spot through the two observation points.
    S_1 = S_0 .* grow_1 .* exp(Df_1);
    S_2 = S_1 .* grow_2 .* exp(Df_2);

    % A path can trigger at most once. Survivors are paths that stayed at
    % or above strike on both observations.
    trig_1  = (S_1 < K);
    trig_2  = (~trig_1) & (S_2 < K);
    survive = ~trig_1 & ~trig_2;

    % PV of the structured coupon leg paid by IB.
    pv_coupon_sim = notional * ( ...
            B_t1_pay * trig_1  * accr_1 * coupon_high + ...
            B_t2_pay * trig_2  * accr_2 * coupon_high + ...
            B_t3_pay * survive * accr_3 * coupon_low );
    pv_coupon = mean(pv_coupon_sim);

    % PV of the floating leg paid by Bank XX. The leg cancels at the trigger
    % year, so the floater identity 1 - B(0, T_end) uses a different T_end
    % per scenario, and so does the BPV that prices the spread.
    sd_dt        = datetime(startDate, 'ConvertFrom', 'datenum');
    t3_pay_unadj = datenum(sd_dt + calyears(3));

    BPV_y1 = BasisPointValueFloating(startDate, ...
                 parameters.firstCouponDate_unadj, d_curve, B_curve);
    BPV_y2 = BasisPointValueFloating(startDate, ...
                 parameters.maturityDate_unadj,    d_curve, B_curve);
    BPV_y3 = BasisPointValueFloating(startDate, t3_pay_unadj, ...
                 d_curve, B_curve);

    pv_A_sim = notional * ( ...
            trig_1  .* ((1 - B_t1_pay) + spread * BPV_y1) + ...
            trig_2  .* ((1 - B_t2_pay) + spread * BPV_y2) + ...
            survive .* ((1 - B_t3_pay) + spread * BPV_y3) );
    pv_A = mean(pv_A_sim);

    upfront = (pv_A - pv_coupon) / notional;
end


function Df = nig_increment(sigma, kappa, eta, alpha, dt, N)
% One NMVM log-return increment over a period of length dt. NIG corresponds
% to alpha = 1/2. Same algebra used for point a, just parametrised on dt
% so we can call it twice for the chained scheme.

    G = random('InverseGaussian', 1, dt / kappa, [N, 1]);
    g = randn(N, 1);

    ln_L_eta = (dt / kappa) * ((1 - alpha) / alpha) * ...
               (1 - (1 + (kappa * eta * sigma^2) / (1 - alpha))^alpha);

    Df = sqrt(dt) * sigma .* sqrt(G) .* g ...
       - (0.5 + eta) * dt * sigma^2 .* G ...
       - ln_L_eta;
end
function upfront = compute_upfront_3y(parameters, sigma, kappa, eta, ...
    N_sim, alpha, coupon_high, coupon_low, spread, notional, ...
    method, sigma_black)

% Mid-market upfront for the 3y version of the certificate (points c/d).
%
% Default (method omitted) uses NIG Monte Carlo and requires (sigma, kappa,
% eta, alpha). Passing method = 'Black' switches to a log-normal Black MC
% and uses sigma_black instead -- the NIG parameters are then ignored.

    if nargin < 11 || isempty(method)
        method = 'NIG';
    end
    useBlack = strcmpi(method, 'Black');
    if useBlack && (nargin < 12 || isempty(sigma_black))
        error('compute_upfront_3y:missingSigmaBlack', ...
              'method = ''Black'' requires sigma_black as 12th argument.');
    end
    rng(42);

    DC_ACT365 = 3;
    DC_30_360 = 6;

    S_0       = parameters.spot;
    div       = parameters.dividend;
    K         = parameters.strike;
    d_curve   = parameters.dates_bootstrap;
    B_curve   = parameters.discounts_bootstrap;
    startDate = parameters.startDate;

    % Payment and reset dates.
    t1_pay   = parameters.firstCouponDate;
    t2_pay   = parameters.maturityDate;
    t3_pay   = business_date_offset(startDate, 'year_offset', 3, ...
                   'convention', 'modified_following');
    t1_reset = parameters.couponResetDate;
    t2_reset = business_date_offset(t2_pay, 'day_offset', -2, ...
                   'convention', 'modified_following');

    % Discount factors.
    B_t1_reset = get_discount_factor_by_zero_rates_linear_interp(startDate, ...
                     t1_reset, d_curve, B_curve);
    B_t2_reset = get_discount_factor_by_zero_rates_linear_interp(startDate, ...
                     t2_reset, d_curve, B_curve);
    B_t1_pay = parameters.B_t1;
    B_t2_pay = parameters.B_t2;
    B_t3_pay = get_discount_factor_by_zero_rates_linear_interp(startDate, ...
                   t3_pay, d_curve, B_curve);

    % Act/365 times driving the dynamics.
    ttm_1 = parameters.ttm_reset;
    ttm_2 = yearfrac(startDate, t2_reset, DC_ACT365);
    dt_1  = ttm_1;
    dt_2  = ttm_2 - ttm_1;

    % 30/360 coupon accruals.
    accr_1 = parameters.accrual_1;
    accr_2 = parameters.accrual_2;
    accr_3 = yearfrac(t2_pay, t3_pay, DC_30_360);

    % Deterministic spot growth between observations.
    grow_1 = exp(-div * dt_1) / B_t1_reset;
    grow_2 = exp(-div * dt_2) * B_t1_reset / B_t2_reset;

    % Two independent log-return increments. Both models are martingale-
    % corrected so that E[exp(Df_i)] = 1 and the deterministic factor
    % grow_i carries the full drift.
    if useBlack
        Df_1 = black_increment(sigma_black, dt_1, N_sim);
        Df_2 = black_increment(sigma_black, dt_2, N_sim);
    else
        Df_1 = nig_increment(sigma, kappa, eta, alpha, dt_1, N_sim);
        Df_2 = nig_increment(sigma, kappa, eta, alpha, dt_2, N_sim);
    end

    % Chain the spot through the two observation points.
    S_1 = S_0 .* grow_1 .* exp(Df_1);
    S_2 = S_1 .* grow_2 .* exp(Df_2);

    trig_1  = (S_1 < K);
    trig_2  = (~trig_1) & (S_2 < K);
    survive = ~trig_1 & ~trig_2;

    fprintf('--- 3Y certificate (%s) ---\n', method);
    fprintf('Q(trigger at T1)              = %.4f\n', mean(trig_1));
    fprintf('Q(trigger at T2 | no T1)      = %.4f\n', sum(trig_2) / sum(~trig_1));
    fprintf('Q(survival to T3)             = %.4f\n', mean(survive));

    % PV of the structured coupon leg paid by IB.
    pv_coupon_sim = notional * ( ...
            B_t1_pay * trig_1  * accr_1 * coupon_high + ...
            B_t2_pay * trig_2  * accr_2 * coupon_high + ...
            B_t3_pay * survive * accr_3 * coupon_low );
    pv_coupon = mean(pv_coupon_sim);

    % PV of the floating leg paid by Bank XX, stops at trigger.
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
% NMVM (NIG when alpha = 1/2) log-return increment over dt, with
% martingale correction so that E[exp(Df)] = 1.

    G = random('InverseGaussian', 1, dt / kappa, [N, 1]);
    g = randn(N, 1);

    ln_L_eta = (dt / kappa) * ((1 - alpha) / alpha) * ...
               (1 - (1 + (kappa * eta * sigma^2) / (1 - alpha))^alpha);

    Df = sqrt(dt) * sigma .* sqrt(G) .* g ...
       - (0.5 + eta) * dt * sigma^2 .* G ...
       - ln_L_eta;
end


function Df = black_increment(sigma, dt, N)
% Black (log-normal) log-return increment over dt, martingale-corrected:
%     Df = -sigma^2/2 * dt + sigma * sqrt(dt) * z,   z ~ N(0,1).

    z  = randn(N, 1);
    Df = -0.5 * sigma^2 * dt + sigma * sqrt(dt) .* z;
end
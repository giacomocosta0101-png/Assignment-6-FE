function dg_price_smile = digital_price_smile(parameters)
% DIGITAL_PRICE_SMILE Calculates the Digital Call price considering the Volatility Smile.
%
% INPUTS:
%   disc_fact    : [Scalar] Discount factor for the maturity
%   volatilities : [Vector] Market volatility surface (specific volatility for each strike)
%   strikes      : [Vector] Array of Strike Prices (K)
%   maturity     : [Datenum] Option expiration date
%   spot_price   : [Scalar] Current spot price (S0)
%   dividends    : [Scalar] Continuous dividend yield
%   discounts    : [Vector] Discount factors extracted from the bootstrap curve
%   startDate    : [Scalar] Start date
%
% OUTPUTS:
%   dg_price_smile : [Vector] Digital prices reflecting market smile
    
    spot_price = parameters.spot;
    rate = parameters.rate; % Extract the risk-free interest rate from parameters
    ttm = parameters.ttm_reset;   % Extract the time to maturity from parameters
    dividends = parameters.dividend; % Extract the continuous dividend yield
    volatilities = parameters.volatilities;
    strikes = parameters.strikes;
    strike = parameters.strike;
    sigma_black = parameters.sigma_black;

    % We calculate Plain Vanilla Call prices for each strike using their 
    % corresponding market volatilities
    call_K = blsprice(spot_price, strike, rate, ttm, sigma_black, dividends);

    % Finite Difference Method:
    % To capture the smile, we perturb the strikes by a small epsilon
    epsilon = 1;
    K_eps = strike+epsilon;

    % We interpolate the volatility surface at the perturbed strikes 
    interp_vol = interp1(strikes, volatilities, K_eps, 'spline', 'extrap');

    % We calculate the price of the Call at perturbed strikes
    call_K_plus_eps = blsprice(spot_price, K_eps, rate, ttm, interp_vol, dividends);

    % We approximate the Digital Price  as [Call(K) - Call(K + epsilon)] / epsilon.
    dg_price_smile = (call_K - call_K_plus_eps) / epsilon;

end
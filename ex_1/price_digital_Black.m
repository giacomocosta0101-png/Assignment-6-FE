function price_digital_put= price_digital_Black(parameters,sigma)
% PRICE_DIGITAL_BLACK Calculates the price of a Digital option using 
% the Black-Scholes model.
% This function assumes a "Flat Volatility" approach by using the ATM volatility 
% for all strikes.
%
% INPUTS:
%   B            : [Scalar] Discount factor for the maturity
%   volatilities : [Vector] Market implied volatilities corresponding to different strikes
%   strikes      : [Vector] Strike prices (K)
%   spot_price   : [Scalar] Current price of the underlying asset (S0) for which we have 
%                           an ATM Forward option
%   dividend     : [Scalar] Continuous annual dividend yield 
%   maturity     : [Datenum] Expiration date of the option
%   discounts    : [Vector] Discount factors extracted from the bootstrap curve
%   startDate    : [Scalar] Start date
%
% OUTPUTS:
%   price_digital_flat_vol : [Vector] Prices of the digital option
%   strike_ATM             : [Scalar] Strike for which we have an ATM
%                                     Forward option
    
    ttm = parameters.ttm_reset;
    F0 = parameters.F0_reset;
    strike = parameters.strike;
    B_reset= parameters.B_reset;

    d2 = (log(F0 ./ strike) - 0.5 * sigma^2 * ttm) / (sigma * sqrt(ttm));

    price_digital_put = B_reset* normcdf(-d2);
    
end
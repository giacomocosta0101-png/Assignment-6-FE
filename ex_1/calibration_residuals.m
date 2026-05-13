function residuals = calibration_residuals(p, x, vol_mkt, t, B, F0,rate,alpha, N, dx)
% INPUTS:
%   p           - [Vector] Current parameters [sigma, kappa, eta]
%   x           - [Vector] Market log-moneyness grid
%   vol_mkt     - [Vector] Market target implied volatilities
%   t           - [Scalar] Time to maturity in years
%   B           - [Scalar] Discount factor for maturity t
%   F0          - [Scalar] Current Forward price
%   alpha       - [Scalar] Tail index (0.5 for NIG, <0.5 for Generalized IG)
%   rate        - [Scalar] Risk-free rate 
%   N           - [Scalar] Number of FFT points
%   dx          - [Scalar] Step size for FFT x-discretization

% OUTPUT:
%   residuals  : [Vector] Difference between model and market prices

    sigma = p(1); % Volatility parameter
    kappa = p(2); % Vol-of-vol parameter 
    eta   = p(3); % Skewness/Drift correction parameter

    % We compute omega_bar that define the boundary of the Laplace exponent domain
    omega_bar = (1 - alpha) / (kappa * sigma^2);
    
    % We check the required condition: eta > -omega_bar 
    if eta <= -omega_bar 
        % Return high residuals if the constraint is not satisfied
        residuals = ones(size(vol_mkt)) * 1e9;
        return;
    end
    
    % Model Pricing using the FFT-based Lewis formula
    prices_model = Price_FFT_NMV(x, sigma, kappa, eta, t, N, dx, B, F0, alpha);

    % We compute the market price by using
    K = F0 .* exp(-x); % strike
    mkt_price = blkprice(F0,K,rate,t,vol_mkt);

    % Assign penalty for NaN results (out of BS bounds) or complex numbers
    prices_model(isnan(prices_model) | ~isreal(prices_model)) = 1e9; 

    % We compute residuals vector for calibration
    residuals = prices_model-mkt_price;
end
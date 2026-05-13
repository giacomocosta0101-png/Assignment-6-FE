function price_FFT_NMV = Price_FFT_NMV(x_grid, sigma, kappa, eta, t, N, dx, B, F0, alpha)
    % PRICE_FFT_NMV computes European option prices using the Fast Fourier 
    % Transform (FFT) based on the Lewis pricing formula for a Normal 
    % Mean-Variance Mixture model.
    %
    % INPUTS:
    %   x_grid          - [Vector] Log-Moneyness values
    %   sigma           - [Scalar] Volatility parameter
    %   kappa           - [Scalar] Vol-of-vol parameter 
    %   eta             - [Scalar] Skewness/Drift correction parameter
    %   t               - [Scalar] Time to maturity in years
    %   N               - [Scalar] Number of FFT points
    %   dx              - [Scalar] Step size for FFT x-discretization
    %   B               - [Scalar] Discount factor for maturity t
    %   F0              - [Scalar] Current Forward price
    %   alpha           - [Scalar] Tail index (0.5 for NIG, <0.5 for Generalized IG)
    %
    % OUTPUT:
    %   price_FFT_NMV   - [Vector] Call option prices for x in x_grid

    % We compute the step size for x-discretization
     dz = (2 * pi) / (N * dx); 
    
    % We define the starting points (z1, x1) 
    z1 = -dz * (N - 1) / 2; 
    x1 = -dx * (N - 1) / 2; 
    
    % We define the grids 
    z = z1 + (0:N-1) * dz;  % Frequency grid
    x = x1 + (0:N-1) * dx;  % Log-Moneyness grid
    
    % Call the model-specific integrand function 
    integrand = integrand_NMV_Lewis(z, t, sigma, kappa, eta, alpha);
    
    % Shift the input because our moneyness grid starts at x1.
    % We also multiply by the integration step 'dz'.
    input_fft = integrand .* exp(-1i .* z .* x1) .* dz;
    
    % Execute the FFT algorithm to compute the discrete sum efficiently
    Y = fft(input_fft);
    
    % Shift the output back because our frequency grid started at z1.
    xk = (x - x1);
    prefactor = exp(-1i .* z1 .* xk);
    integral_vals = prefactor .* Y;

    % Apply the Lewis formula and take the real part to remove negligible 
    % imaginary numerical noise.
    price_FFT = F0 .* B .* (1 - (exp(-x / 2)) / (2 * pi) .* real(integral_vals));
    
    % FFT results can be numerically unstable in the far tails.
    % We filter the x-grid (e.g., between -1 and 1) to remove noise before interpolation.
    idx_valid = find(x > -1 & x < 1);
    x_clean = x(idx_valid);
    price_FFT_clean = price_FFT(idx_valid);

    % Interpolate FFT results to match the requested log-moneyness grid
    % (x_grid)
    price_FFT_NMV = interp1(x_clean, price_FFT_clean, x_grid, 'spline');
end
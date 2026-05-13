function integrand_NMV_Lewis = integrand_NMV_Lewis (z, t, sigma, k, eta, alpha)

    % INTEGRAND_NMV_LEWIS computes the core of the Lewis integral 
    % for the Normal Tempered Stable Mean-Variance Mixture model.
    %
    % INPUT:
    %   z                       : [Vector] Integration variable 
    %   t                       : [Scalar] Time to maturity in years 
    %   sigma                   : [Scalar] Volatility parameter 
    %   k                       : [Scalar] Vol-of-vol parameter 
    %   eta                     : [Scalar] Skewness parameter
    %   alpha                   : [Scalar] Tempered stable parameter, must be in (0,1)
    %
    % OUTPUT:
    %   integrand_NMV_Lewis      : [Vector] Complex integrand values for the Lewis formula

        
    % Generalized Laplace Exponent
    ln_Lap = @(w) (t / k) * ((1 - alpha) / alpha) .* (1 - (1 + (w .* k .* sigma^2) ./ (1 - alpha)).^alpha);
    
    % Martingale drift correction component (evaluated at eta)
    ln_Lap_eta = ln_Lap(eta);
    
    % Frequency shift for the Lewis formula
    u = -z - 0.5i;
    
    % Argument for the Laplace transform inside the Characteristic Function
    argum_Lap = (u.^2 + 1i * (1 + 2 * eta) .* u) / 2;
    
    % Assemble the Characteristic Function phi(u)
    ln_Lap_argum = ln_Lap(argum_Lap);
    phi = exp(-1i * u .* ln_Lap_eta + ln_Lap_argum);
    
    % Lewis Integrand Core
    integrand_NMV_Lewis = phi ./ (z.^2 + 0.25);
    
end
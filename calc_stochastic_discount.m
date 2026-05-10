function [D_u, D_m, D_d] = calc_stochastic_discount(x_nodes, B_t_t1, dx, mu_hat, a, sigma, dt)
% CALC_STOCHASTIC_DISCOUNT Computes the exact stochastic discount factors 
% for the Hull-White trinomial tree between t_i and t_{i+1}.
%
% INPUTS:
%   x_nodes     : Column vector of spatial grid nodes at current step (x_i)
%   B_t_t1      : Deterministic forward discount factor B(t_i, t_{i+1})
%   dx          : Spatial step size of the tree (Delta x)
%   mu_hat      : Exact drift coefficient (1 - exp(-a*dt))
%
% OUTPUTS:
%   D_u, D_m, D_d : Column vectors of stochastic discounts for UP, MID, DOWN

sigma_hat = sigma * sqrt((1 - exp(-2 * a * dt)) / (2 * a));

term1 = dt;
term2 = 2 * (1 - exp(-a * dt)) / a;
term3 = (1 - exp(-2 * a * dt)) / (2 * a);
sigma_star = (sigma / a) * sqrt(term1 - term2 + term3);

% Pre-computed terms for the exponential formula 
term_const = -0.5 * (sigma_star^2);
ratio_sigma = sigma_star / sigma_hat;

N_nodes = length(x_nodes);

% --- 1. JUMPS DEFINITION (Delta x_{i+1}) ---
% We define the physical jump size the process makes on the grid
dx_u = repmat(dx, N_nodes, 1);    % UP branch jumps by +dx
dx_m = zeros(N_nodes, 1);         % MID branch stays on same level (jump = 0)
dx_d = repmat(-dx, N_nodes, 1);   % DOWN branch jumps by -dx

% Boundary CASE B: Bottom of the tree (Row 1)
dx_d(1) = 0;        % Lowest branch stays at l
dx_m(1) = dx;       % Middle branch jumps up 1 level
dx_u(1) = 2 * dx;   % Highest branch jumps up 2 levels

% Boundary CASE C: Top of the tree (Last row)
dx_u(end) = 0;      % Highest branch stays at l
dx_m(end) = -dx;    % Middle branch jumps down 1 level
dx_d(end) = -2 * dx;% Lowest branch jumps down 2 levels


% --- 2. EXACT SLIDE FORMULA APPLICATION ---
% Compute the deterministic drift term for all nodes: [\hat{\mu} * x_i]
drift_term = mu_hat .* x_nodes;

% Apply the exact formula from the slide: 
% exp{ -0.5*(sigma_star)^2 - (sigma_star/sigma_hat) * [Delta_x + mu_hat*x_i] }
% Multiplied by the forward discount factor B(t_i, t_{i+1})

D_u = B_t_t1 * exp(term_const - ratio_sigma * (dx_u + drift_term));
D_m = B_t_t1 * exp(term_const - ratio_sigma * (dx_m + drift_term));
D_d = B_t_t1 * exp(term_const - ratio_sigma * (dx_d + drift_term));

end
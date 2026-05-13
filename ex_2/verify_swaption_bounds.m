function validation_table = verify_swaption_bounds(precision_levels, prices, ...
    LB, UB_eur, UB_cap)
% VERIFY_SWAPTION_BOUNDS Validation table for the Bermudan tree price vs analytical bounds.
%
% Reports three bounds for each precision level:
%   LB      : tight lower bound (max European on residual swap)
%   UB_cap  : tight upper bound (sum of single-period caplets) -- main check
%   UB_eur  : loose upper bound (sum of Europeans on residual swap) -- for reference
%
% The OK/FAIL check uses the tight bounds [LB, UB_cap]. If the tree price
% violates this interval, the tree implementation is likely inconsistent.
%
% INPUTS:
%   precision_levels : [vector] Steps-per-year used in the convergence test.
%   prices           : [vector] Bermudan prices from the tree.
%   LB               : [scalar] Lower bound (Jamshidian max European).
%   UB_eur           : [scalar] Loose upper bound (Jamshidian sum of Europeans).
%   UB_cap           : [scalar] Tight upper bound (sum of single-period caplets).
%
% OUTPUT:
%   validation_table : [table] Summary table with all three bounds and check.

    n_levels = length(precision_levels);

    Lower_Bound      = repmat(LB,     n_levels, 1);
    Upper_Cap        = repmat(UB_cap, n_levels, 1);
    Upper_European   = repmat(UB_eur, n_levels, 1);

    % Consistency check on the TIGHT interval [LB, UB_cap]
    is_ok = (prices >= LB - 1e-9) & (prices <= UB_cap + 1e-9);

    Validation = repmat("FAIL - Out of Bounds", n_levels, 1);
    Validation(is_ok) = "OK - Consistent";

    validation_table = table(precision_levels(:), prices(:), ...
        Lower_Bound, Upper_Cap, Upper_European, Validation, ...
        'VariableNames', {'Steps_Per_Year', 'Tree_Price', ...
                          'Lower_Bound', 'Upper_Bound_Caplet', ...
                          'Upper_Bound_European', 'Check'});
end

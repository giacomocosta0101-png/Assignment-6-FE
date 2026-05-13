function [time_grid, grid_dates, dt_vector, idx_exercise] = build_hw_time_grid(...
    startDate, maturityDate, stepsPerYear, exerciseDates)
% BUILD_HW_TIME_GRID Piecewise-uniform time grid for HW Bermudan/European/ZCB pricing.
%
% APPROACH:
% Internal tree nodes are NOT business-day-adjusted. Only contractual dates
% (start, exercise dates, maturity) are guaranteed to be tree nodes. Within
% each segment between two consecutive contractual dates, the sub-grid is
% strictly uniform with approximately `stepsPerYear` steps per year.
%
% This ensures that dt is EXACTLY constant within each segment, which is
% what the Hull-White closed-form formulas (probabilities, stochastic
% discounts, sigma_hat, sigma_star) assume. dt is allowed to vary across
% segments because each backward-induction step uses its own dt_vector(i).
%
% INPUTS:
%   startDate     : valuation/start date (datetime or datenum)
%   maturityDate  : final maturity (datetime or datenum)
%   stepsPerYear  : approximate number of tree steps per year (scalar)
%   exerciseDates : (optional) contractual exercise dates.
%                   Omit or pass [] for instruments without early exercise (ZCB).
%
% OUTPUTS:
%   time_grid    : [N x 1] times in years (Act/365), strictly increasing
%   grid_dates   : [N x 1] datenums = startDate + time_grid * 365 (NOT BD-adjusted)
%   dt_vector    : [(N-1) x 1] intervals between consecutive nodes
%   idx_exercise : [n_ex x 1] indices in time_grid corresponding to exerciseDates
%                  (empty if no exercise dates given)

    % 1. UNIFY INPUT TYPES
    startDate    = datenum(startDate);
    maturityDate = datenum(maturityDate);
    
    if nargin < 4 || isempty(exerciseDates)
        T_ex_y = [];
    else
        T_ex_y = yearfrac(startDate, datenum(exerciseDates), 3);
        T_ex_y = T_ex_y(:);                             % force column
    end
    
    % 2. TOTAL HORIZON
    T_max_y = yearfrac(startDate, maturityDate, 3);
    
    % 3. ANCHOR TIMES: start, exercise times, maturity (deduplicated, sorted)
    anchors = unique([0; T_ex_y; T_max_y]);
    
    % 4. PIECEWISE-UNIFORM SUB-GRIDS
    % Within each segment [anchor_i, anchor_{i+1}] build a uniform partition
    % with ~stepsPerYear steps per year (at least 1 step per segment).
    time_grid = anchors(1);                             % start with t = 0
    for i = 1:length(anchors) - 1
        seg_len  = anchors(i+1) - anchors(i);
        n_sub    = max(1, round(seg_len * stepsPerYear));   % # of sub-steps
        sub_grid = linspace(anchors(i), anchors(i+1), n_sub + 1)';
        time_grid = [time_grid; sub_grid(2:end)];       % append, skip duplicated boundary
    end
    
    % 5. dt VECTOR
    dt_vector = diff(time_grid);
    
    % 6. LOCATE EXERCISE DATES IN THE FINAL GRID
    if isempty(T_ex_y)
        idx_exercise = [];
    else
        tol = 1e-10;
        idx_exercise = arrayfun(@(t) find(abs(time_grid - t) < tol, 1), T_ex_y);
    end
    
    % 7. CORRESPONDING DATENUMS (linear Act/365 conversion)
    % These are valid floating-point datenums: yearfrac and curve interpolation
    % handle them correctly. They are NOT business-day-adjusted by design.
    grid_dates = startDate + time_grid * 365;

end

function [time_grid_y, grid_dates, dt_vector, idx_exercise] = build_hw_time_grid(...
    startDate, maturityDate, exerciseDates, stepsPerYear)
% BUILD_HW_TIME_GRID Constructs a time grid for Hull-White Bermudan Swaption pricing.
%
% This function creates a uniform time grid, adjusts nodes to business days,
% and injects exact exercise dates into the grid ensuring no duplicates.
%
% INPUTS:
%   startDate      : [datetime] Valuation/Start date.
%   maturityDate   : [datetime] Final maturity of the underlying swap.
%   exerciseDates  : [datetime vector] Dates when the option can be exercised.
%   stepsPerYear   : [scalar] Desired number of steps per year.
%
% OUTPUTS:
%   time_grid_y    : [vector] Time steps in years (Act/365).
%   grid_dates     : [datetime vector] Corresponding business-adjusted dates.
%   dt_vector      : [vector] Time intervals (dt) between nodes.
%   idx_exercise   : [vector] Indices in grid_dates corresponding to exerciseDates.

    % 1. Initial Uniform Grid Construction
    T_max_y = yearfrac(startDate, maturityDate, 3);
    num_steps = max(1, round(T_max_y * stepsPerYear));
    dt_init = T_max_y / num_steps;
    
    % Raw time grid
    time_grid_y = (0:num_steps)' * dt_init;
    
    % 2. Convert to Dates and Apply Business Day Adjustment
    grid_dates_raw = startDate + years(time_grid_y);
    grid_dates_num = arrayfun(@(d) business_date_offset(d), grid_dates_raw);
    grid_dates = datetime(grid_dates_num, 'ConvertFrom', 'datenum')';

    % 3. Inject Exercise Dates
    T_ex_y = yearfrac(startDate, exerciseDates, 3);
    
    % Combine initial grid with exercise dates
    combined_T = [time_grid_y; T_ex_y(:)];
    combined_D = [grid_dates(:); exerciseDates(:)];
    
    % Sort by time
    [combined_T, sortIdx] = sort(combined_T);
    combined_D = combined_D(sortIdx);
    
    % Vectorized De-duplication (Snapping)
    % If two nodes are closer than 1 hour, we keep the one coming from exerciseDates
    tol = 1 / (365 * 24); 
    diff_T = diff(combined_T);
    to_remove = false(size(combined_T));
    
    % Find nodes too close to each other
    close_idx = find(diff_T < tol);
    
    % In case of proximity, we mark the "non-exercise" node for removal
    % We prioritize exerciseDates which are at the end of the original combined arrays
    to_remove(close_idx) = true; 
    
    time_grid_y = combined_T(~to_remove);
    grid_dates = combined_D(~to_remove);
    
    % 4. Final Outputs
    dt_vector = diff(time_grid_y);
    
    % Map exercise dates to the final grid 
    % We use min(abs) for each exercise date to find its position in the refined grid
    idx_exercise = arrayfun(@(ex) find(abs(time_grid_y - ex) < tol, 1), T_ex_y);

end
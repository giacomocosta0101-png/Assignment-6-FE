function print_case_study_1(notional, N_sim, ...
    X_nig_2y, X_black_smile_2y, X_nig_3y, ...
    X_black_flat_2y, X_black_flat_3y)
% Pretty-print upfronts for the certificate swap (points a, b, d, e).

    line = repmat('-', 1, 64);
    row  = @(label, X) fprintf(' %-22s %10.4f %14.2f %14.0f\n', ...
                label, X*100, X*1e4, X*notional);
    diff_row = @(label, dX) fprintf('  %-3s %+8.4f %%   ( %+7.2f bps, %+12.0f EUR )\n', ...
                label, dX*100, dX*1e4, dX*notional);

    fprintf('\n%s\n', line);
    fprintf(' Case Study 1 -- Certificate swap upfronts\n');
    fprintf(' Notional = %.0f Mln EUR,  N_sim = %.0e\n', notional/1e6, N_sim);
    fprintf('%s\n', line);
    fprintf(' %-22s %10s %14s %14s\n', 'Model', 'X [%]', 'X [bps]', 'EUR');
    fprintf('%s\n', line);

    fprintf(' 2-year contract\n');
    row('  NIG  (point a)',     X_nig_2y);
    row('  Black smile (b)',    X_black_smile_2y);
    row('  Black flat   (e)',   X_black_flat_2y);

    fprintf('\n 3-year contract\n');
    row('  NIG  (point d)',     X_nig_3y);
    row('  Black flat   (e)',   X_black_flat_3y);

    fprintf('%s\n', line);
    fprintf(' Model error (Black flat - NIG)\n');
    diff_row('2Y:', X_black_flat_2y - X_nig_2y);
    diff_row('3Y:', X_black_flat_3y - X_nig_3y);
    fprintf('%s\n\n', line);
end
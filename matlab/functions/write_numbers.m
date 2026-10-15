function write_numbers(cfg, name, rows)
%WRITE_NUMBERS Save numbers quoted in the paper text to output/tables/numbers_<name>.csv.
%   rows: cell array, one row per number: {item, value, lo, hi, note};
%   lo and hi are NaN when the number has no interval.

tbl = cell2table(rows, 'VariableNames', {'item', 'value', 'lo', 'hi', 'note'});
f = fullfile(cfg.table_dir, ['numbers_' name '.csv']);
writetable(tbl, f);
fprintf('Saved: %s\n', f);
end

function s = format_summary(v, scale, fmt, interval)
%FORMAT_SUMMARY Text "median [lo, hi]" of draws v (scaled by scale).
%   interval: 'ci95' (2.5-97.5%, default) or 'iqr' (25-75%)
%   fmt     : sprintf format with three numeric fields

if nargin < 4; interval = 'ci95'; end
switch interval
    case 'ci95'; p = [0.025 0.975];
    case 'iqr';  p = [0.25 0.75];
    otherwise;   error('interval must be ci95 or iqr');
end
s = sprintf(fmt, scale * median(v), scale * quantile(v, p(1)), scale * quantile(v, p(2)));
end

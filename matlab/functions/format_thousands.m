function s = format_thousands(x)
%FORMAT_THOUSANDS Round to an integer and insert thousands separators.

s = sprintf('%.0f', round(x));
neg = startsWith(s, '-');
if neg; s = s(2:end); end
first = mod(numel(s), 3);
if first == 0; first = 3; end
parts = {s(1:first)};
for i = first+1 : 3 : numel(s)
    parts{end+1} = s(i:i+2); %#ok<AGROW>
end
s = strjoin(parts, ',');
if neg; s = ['-' s]; end
end

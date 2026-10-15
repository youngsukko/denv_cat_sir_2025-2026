function shp = load_brazil_shapes(cfg)
%LOAD_BRAZIL_SHAPES State polygons (IBGE BR_UF_2023) with plotting helpers.
%   shp(i).code  state abbreviation
%   shp(i).si    index into cfg.state_codes
%   shp(i).parts cell array of [x y] rings with at least three points
%   shp(i).cx/cy label position (centroid of all vertices as one polyshape)

raw = shaperead(cfg.file_shape);
shp = struct('code', {}, 'si', {}, 'parts', {}, 'cx', {}, 'cy', {});
for i = 1:numel(raw)
    x = raw(i).X(:);
    y = raw(i).Y(:);
    brk = [0; find(isnan(x)); numel(x) + 1];
    parts = {};
    for k = 1:numel(brk) - 1
        seg = brk(k)+1 : brk(k+1)-1;
        if numel(seg) >= 3; parts{end+1} = [x(seg) y(seg)]; end %#ok<AGROW>
    end
    valid = ~isnan(x) & ~isnan(y);
    [cx, cy] = centroid(polyshape(x(valid), y(valid)));

    shp(i).code  = raw(i).SIGLA_UF;
    shp(i).si    = find(strcmp(cfg.state_codes, raw(i).SIGLA_UF));
    shp(i).parts = parts;
    shp(i).cx    = cx;
    shp(i).cy    = cy;
end
end
